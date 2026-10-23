module MetricCodifferentialTests

using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Fields
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.TensorValues
using FillArrays
using Test
import LinearAlgebra

# The codifferential on a chart model: `δ` reads the metric of the triangulation
# and the Christoffel symbols its derivatives give. The `⋆d⋆` spelling of the
# same operator and the adjointness to `𝑑` are what check it.

# ── Chart models of given metric ──────────────────────────────────────────────

struct ChartModel{Dc,Dp,M,A} <: DiscreteModel{Dc,Dp}
  model::M
  cell_metric::A
  function ChartModel(model::DiscreteModel{Dc,Dp},cell_metric) where {Dc,Dp}
    new{Dc,Dp,typeof(model),typeof(cell_metric)}(model,cell_metric)
  end
end

Gridap.Geometry.get_grid(m::ChartModel)          = get_grid(m.model)
Gridap.Geometry.get_grid_topology(m::ChartModel) = get_grid_topology(m.model)
Gridap.Geometry.get_face_labeling(m::ChartModel) = get_face_labeling(m.model)
Gridap.Geometry.get_cell_metric(m::ChartModel)   = m.cell_metric

chart_model(model,g) = ChartModel(model,Fill(GenericField(g),num_cells(model)))

# Covariant as a metric must be, non constant so that the Christoffel symbols do
# not vanish, and diagonally dominant hence positive definite, which is the
# metric index 0 that the sign of `±⋆d⋆` below assumes.
g2(x) = SymTensorValue{Cov(2)}(1.5 + 0.4x[1]^2, 0.25x[2], 2.0 + 0.3x[1]*x[2])
g3(x) = SymTensorValue{Cov(3)}(1.5 + 0.4x[1]^2, 0.25x[2], 0.1x[3],
                               2.0 + 0.3x[1]*x[2], 0.2x[1], 1.75 + 0.5x[3]^2)

model2 = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(2,2))
model3 = CartesianDiscreteModel((0.0,1.0,0.0,1.0,0.0,1.0),(2,2,2))
cmodel2 = chart_model(model2, g2)
cmodel3 = chart_model(model3, g3)

# A K-form whose every component varies in every coordinate
mk(K,D) = x -> ExteriorFormValue{K,Cov(D)}(
  ntuple(i -> 0.3 + 0.7i*x[1] + 0.4x[2]^2 - 0.2i*x[mod1(i,D)]*x[D], binomial(D,K)))

form_field(K,D,trian) = CellField(mk(K,D), trian, PhysicalDomain())

# The largest difference between `a` and `sgn*b`, two cell arrays of form values,
# over every cell and every point.
function max_diff(a, b, sgn=1)
  ca, cb = array_cache(a), array_cache(b)
  m = 0.0
  for cell in eachindex(a,b)
    va, vb = getindex!(ca,a,cell), getindex!(cb,b,cell)
    for q in eachindex(va,vb)
      m = max(m, maximum(abs, Tuple(va[q]) .- sgn.*Tuple(vb[q]); init=0.0))
    end
  end
  m
end

# ── `δ` of a chart is `±⋆d⋆` of the same metric ───────────────────────────────

# `δα = (-1)^{D(K+1)+1} ⋆d⋆α` for a positive definite metric. The two sides
# share no kernel: `δ` contracts the Christoffel symbols of the metric, while
# `⋆d⋆` composes two metric Hodge stars with the metric free `𝑑`.
@testset "δ is ±⋆d⋆ on a chart of dimension $D" for (D,cmodel) in
    ((2,cmodel2), (3,cmodel3))
  Ω = Triangulation(cmodel)
  pts = get_cell_points(Ω)
  for K in 1:D
    ω = form_field(K,D,Ω)
    sgn = isodd(D*(K+1)+1) ? -1 : 1
    sdsω = hodge_star(𝑑(hodge_star(ω)))(pts)
    @test max_diff(δ(ω)(pts), sdsω, sgn) < 1.0e-12
    # the sign is the definition's, so the opposite one is nowhere near
    @test max_diff(δ(ω)(pts), sdsω, -sgn) > 1.0e-3
  end
end

# ── The metric the operator reads ─────────────────────────────────────────────

# A chart whose metric is the identity has vanishing symbols, so it answers what
# the euclidean kernel does on the stock model
let imodel = chart_model(model2, x -> SymTensorValue{Cov(2)}(1.0,0.0,1.0))
  Ωi, Ωs = Triangulation(imodel), Triangulation(model2)
  for K in 1:2
    @test max_diff(δ(form_field(K,2,Ωi))(get_cell_points(Ωi)),
                   δ(form_field(K,2,Ωs))(get_cell_points(Ωs))) < 1.0e-14
  end
end

# `LinearAlgebra.I` asks for the euclidean codifferential on a chart model all
# the same, and that is not the chart's
let Ω = Triangulation(cmodel2), Ωs = Triangulation(model2)
  ω, ωs = form_field(1,2,Ω), form_field(1,2,Ωs)
  pts, pts_s = get_cell_points(Ω), get_cell_points(Ωs)
  @test max_diff(δ(ω,LinearAlgebra.I)(pts), δ(ωs)(pts_s)) < 1.0e-14
  @test max_diff(δ(ω)(pts), δ(ω,LinearAlgebra.I)(pts)) > 1.0e-3
end

# A metric given explicitly is the one used, in place of the triangulation's
let Ω = Triangulation(model2), Ωc = Triangulation(cmodel2)
  g = CellField(g2, Ω, PhysicalDomain())
  @test max_diff(δ(form_field(1,2,Ω),g)(get_cell_points(Ω)),
                 δ(form_field(1,2,Ωc))(get_cell_points(Ωc))) < 1.0e-12
end

# ── `δ` is the adjoint of `𝑑` for the metric of the chart ─────────────────────

# `∫(𝑑α|ω) = ∫(α|δω)` for an `α` vanishing on the boundary, both integrals
# taking the volume element of the metric: the property that defines `δ`, and
# the one the third Christoffel term is needed for at `K ≥ 2`.
bump(x) = x[1]*(1-x[1])*x[2]*(1-x[2])

let cmodel = chart_model(CartesianDiscreteModel((0.0,1.0,0.0,1.0),(4,4)), g2)
  Ω = Triangulation(cmodel)
  dΩ = Measure(Ω,8)

  α0 = CellField(x -> ExteriorFormValue{0,Cov(2)}(bump(x)), Ω, PhysicalDomain())
  ω1 = form_field(1,2,Ω)
  @test sum(∫( alt_inner(𝑑(α0), ω1) )dΩ) ≈ sum(∫( alt_inner(α0, δ(ω1)) )dΩ)

  α1 = CellField(x -> ExteriorFormValue{1,Cov(2)}(bump(x)*(1.0 + x[1]),
                                                  bump(x)*(0.5 - x[2])),
                 Ω, PhysicalDomain())
  ω2 = form_field(2,2,Ω)
  @test sum(∫( alt_inner(𝑑(α1), ω2) )dΩ) ≈ sum(∫( alt_inner(α1, δ(ω2)) )dΩ)
end

# ── The derivatives of the metric reach the operators ─────────────────────────

# The inverse metric and the signed volume density are built from the metric, so
# the rules of `inv` and of `sgn(det g)√|det g|` are what lets them be
# differentiated, and with them the metric `⋆` of an expression `𝑑` is taken of.
let Ω = Triangulation(cmodel2), pts = get_cell_points(Ω)
  g, g_inv, signed_sqrt_g = get_metric_data(Ω)
  ∂g, ∂g_inv, ∂signed_sqrt_g = map(partial_derivatives, (g, g_inv, signed_sqrt_g))
  gv, giv, ∂gv, ∂giv, ∂sv = map(
    a -> first(first(a(pts))), (g, g_inv, ∂g, ∂g_inv, ∂signed_sqrt_g))

  @test ∂gv isa ProductTensorValue{SymNTensorValue{Tuple{Cov(2)}},SymTensorValue{Cov(2)}}
  @test ∂giv isa ProductTensorValue{SymNTensorValue{Tuple{Cov(2)}},SymTensorValue{Con(2)}}
  @test ∂sv isa ProductTensorValue{SymNTensorValue{Tuple{Cov(2)}},ScalarValue}

  # `∂(g⁻¹) = -g⁻¹ (∂g) g⁻¹` and `∂√|det g| = √|det g| g^{ab} (∂g_{ab}) / 2`
  for c in 1:2, a in 1:2, b in 1:2
    @test ∂giv[c,a,b] ≈ -sum(giv[a,i]*∂gv[c,i,j]*giv[j,b] for i in 1:2, j in 1:2)
  end
  for c in 1:2
    @test ∂sv[c] ≈
      sqrt(abs(det(gv))) * sum(giv[a,b]*∂gv[c,a,b] for a in 1:2, b in 1:2) / 2
  end
end

# ── A facet of a chart model carries no metric of its own ─────────────────────

# `δ` needs the derivatives of the metric induced on the facet, which the model
# does not hold, so it fails there rather than take the cell metric
let Γ = BoundaryTriangulation(cmodel2)
  @test_throws ErrorException δ(form_field(1,2,Γ))(get_cell_points(Γ))
end

# ── `δ` of a form valued FE basis assembles ───────────────────────────────────

# The basis reaches the kernel through the operation tree, so `δ` of it enters a
# bilinear form like any other operator of the calculus, and the matrix it
# assembles is the chart's rather than the euclidean one.
let tri = simplexify(CartesianDiscreteModel((0,1,0,1),(2,2)))
  ctri = chart_model(tri, g2)
  reffe = ReferenceFE(TRI,:P⁻,1,1,Float64)
  spaces(model) = (V = FESpace(model,reffe); (TrialFESpace(V), V))
  Uc, Vc = spaces(ctri)
  U , V  = spaces(tri)
  dΩc, dΩ = Measure(Triangulation(ctri),4), Measure(Triangulation(tri),4)

  Mc = assemble_matrix((u,v) -> ∫( alt_inner(δ(u),δ(v)) )dΩc, Uc, Vc)
  M  = assemble_matrix((u,v) -> ∫( alt_inner(δ(u),δ(v)) )dΩ , U , V )
  @test Mc ≈ Mc'
  @test !iszero(Mc)
  @test !(Mc ≈ M)
end

end # module
