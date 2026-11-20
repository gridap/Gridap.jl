module CovariantDerivativeTests

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

# `cov∇` on a chart model reads the metric of the triangulation and the
# Christoffel symbols its derivatives give, like `δ` does. What is checked here
# is that plumbing: that the metric reaching the kernel is the triangulation's,
# that `I` overrides it, and that an FE basis reaches it through the operation
# tree. The identities of the connection itself are algebraic in `∂a` and `Γ`,
# and are checked of arbitrary ones in `TensorValuesTests`.

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

g2(x) = SymTensorValue{Cov(2)}(1.5 + 0.4x[1]^2, 0.25x[2], 2.0 + 0.3x[1]*x[2])

model2  = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(2,2))
cmodel2 = chart_model(model2, g2)
Ωc, Ωs  = Triangulation(cmodel2), Triangulation(model2)
ptsc, ptss = get_cell_points(Ωc), get_cell_points(Ωs)

mk(K,D) = x -> ExteriorFormValue{K,Cov(D)}(
  ntuple(i -> 0.3 + 0.7i*x[1] + 0.4x[2]^2 - 0.2i*x[mod1(i,D)]*x[D], binomial(D,K)))

vec_field(D) = x -> ContraVectorValue{Con(D)}(
  ntuple(i -> 0.2 + 0.5i*x[1] - 0.3x[mod1(i+1,D)]^2, D))

form_field(K,D,trian) = CellField(mk(K,D), trian, PhysicalDomain())

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

max_abs(a) = max_diff(a, a, 0)

# ── The metric the operator reads is the triangulation's ──────────────────────

# `∇g = 0` holds of the Christoffel symbols of `g` only, so it vanishes here
# exactly when the metric fed to the kernel is the chart's own. Both variance
# branches are annihilated.
let gfield = CellField(g2, Ωc, PhysicalDomain())
  @test max_abs(cov∇(gfield)(ptsc)) < 1.0e-12
  @test max_abs(cov∇(CellField(x -> inv(g2(x)), Ωc, PhysicalDomain()))(ptsc)) < 1.0e-12
  # and the Christoffel terms are not negligible to begin with
  @test max_abs(∂(gfield)(ptsc)) > 1.0e-2
end

# A chart whose metric is the identity has vanishing symbols, so `∇` answers
# what `∂` does on the stock model
let Ωi = Triangulation(chart_model(model2, x -> SymTensorValue{Cov(2)}(1.0,0.0,1.0)))
  @test max_diff(cov∇(form_field(1,2,Ωi))(get_cell_points(Ωi)),
                 ∂(form_field(1,2,Ωs))(ptss)) < 1.0e-14
end

# The stock model carries no metric, so `∇` is `∂` there without being asked
@test max_diff(cov∇(form_field(1,2,Ωs))(ptss), ∂(form_field(1,2,Ωs))(ptss)) < 1.0e-14

# `LinearAlgebra.I` asks for the euclidean covariant derivative on a chart model
# all the same, in both arities, and that is not the chart's
let ω = form_field(1,2,Ωc), u = CellField(vec_field(2), Ωc, PhysicalDomain())
  @test max_diff(cov∇(ω,LinearAlgebra.I)(ptsc), ∂(form_field(1,2,Ωs))(ptss)) < 1.0e-14
  @test max_diff(cov∇(ω)(ptsc), cov∇(ω,LinearAlgebra.I)(ptsc)) > 1.0e-3
  @test max_diff(cov∇(u,ω,LinearAlgebra.I)(ptsc), cov∇(u,ω)(ptsc)) > 1.0e-3
end

# ── `∇` of a form valued FE basis assembles ───────────────────────────────────

# The basis reaches the kernel through the operation tree, so `∇_v` of it enters
# a bilinear form like any other operator of the calculus.
let tri = simplexify(CartesianDiscreteModel((0,1,0,1),(2,2)))
  ctri = chart_model(tri, g2)
  V = FESpace(ctri, ReferenceFE(TRI,:P⁻,1,1,Float64))
  Ω = Triangulation(ctri)
  dΩ = Measure(Ω,4)
  w = CellField(vec_field(2), Ω, PhysicalDomain())

  M = assemble_matrix((u,v) -> ∫( alt_inner(cov∇(w,u),cov∇(w,v)) )dΩ, TrialFESpace(V), V)
  @test M ≈ M'
  @test !iszero(M)
end

# `cov∇` is the alias of `covariant_derivative`
@test cov∇ === covariant_derivative

end # module
