module MetricFormAssemblyTests

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

# The metric operators and the measure convention through the FE path: a 𝓟ᵣΛᵏ
# space on a simplicial chart model. `∫( a ∧ ⋆b )dΩ` and `∫( alt_inner(a,b) )dΩ`
# are the same bilinear form read from either end, and each is the metric weak
# form written by hand.

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

# Covariant as a metric must be, and constant so that the weak forms written by
# hand below have a closed form.
const g2 = SymTensorValue{Cov(2)}(2.0,0.5,1.0)
const g3 = SymTensorValue{Cov(3)}(2.0,0.5,0.0,1.0,0.25,1.5)

# `g⁻¹` spelled without variance, which a plain `VectorValue` contracts with as
# it stands: a hand written integrand raises form indices itself rather than
# through `sharp`.
unvarianced_inv(g) = TensorValue(inv(get_array(g)))

# The moment based elements need every cell to list its vertices in increasing
# global id order, which `simplexify` of a Cartesian model gives.
sorted_tri = simplexify(CartesianDiscreteModel((0,1,0,1),(2,2)))
sorted_tet = simplexify(CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2)))
ctri = chart_model(sorted_tri, x -> g2)
ctet = chart_model(sorted_tet, x -> g3)

Ωc = Triangulation(ctri); Ω = Triangulation(sorted_tri)
fe(p,F,r,k) = ReferenceFE(p,F,r,k,Float64)
spaces(model,reffe) = (V = FESpace(model,reffe); (TrialFESpace(V), V))

# ── The two spellings of a metric bilinear form assemble the same matrix ──────

# `a ∧ ⋆b = (a|b) √|det g| dx¹∧…∧dx^D`: the wedge is a top form, integrated
# against the euclidean weights because its `⋆` carries the volume factor, and
# `alt_inner` a scalar, integrated against the metric volume element. Both read
# the metric of the model, so the two assemble the same matrix.
@testset "∧⋆ and alt_inner agree on $lbl" for (lbl,cmodel,p,D,cases) in (
    ("TRI", ctri, TRI, 2, ((:P,2,0),(:P⁻,1,1),(:P,2,1),(:P⁻,1,2),(:P,2,2))),
    ("TET", ctet, TET, 3, ((:P⁻,1,0),(:P⁻,1,1),(:P⁻,1,2),(:P⁻,1,3))),
  )
  for (F,r,k) in cases
    U, V = spaces(cmodel, fe(p,F,r,k))
    dΩ = Measure(Triangulation(cmodel), 2r+2)

    M_alt = assemble_matrix((u,v) -> ∫( alt_inner(u,v) )dΩ, U, V)
    @test M_alt ≈ assemble_matrix((u,v) -> ∫( u ∧ hodge_star(v) )dΩ, U, V)

    # The tensor inner product counts each independent component `K!` times.
    @test assemble_matrix((u,v) -> ∫( inner(u,v) )dΩ, U, V) ≈ factorial(k) * M_alt

    k == D && continue
    K_alt = assemble_matrix((u,v) -> ∫( alt_inner(𝑑(u),𝑑(v)) )dΩ, U, V)
    @test K_alt ≈ assemble_matrix((u,v) -> ∫( 𝑑(u) ∧ hodge_star(𝑑(v)) )dΩ, U, V)
  end
end

# `𝑑` of a top form is zero, and the `(D+1)`-forms it is valued in have no `⋆`,
# so that one bilinear form has the `alt_inner` spelling alone.
let dΩ = Measure(Ωc,4), (U,V) = spaces(ctri, fe(TRI,:P⁻,1,2))
  @test iszero(assemble_matrix((u,v) -> ∫( alt_inner(𝑑(u),𝑑(v)) )dΩ, U, V))
  @test_throws ErrorException assemble_matrix((u,v) -> ∫( 𝑑(u) ∧ hodge_star(𝑑(v)) )dΩ, U, V)
end

# ── The assembled form is the metric weak form written by hand ────────────────

# The same space on the stock model carries the same DOFs and the euclidean
# weights, so a hand written integrand there supplies `√|det g|` itself.
@testset "the hand written metric weak form of $lbl" for
    (lbl,cmodel,model,p,D,g,cases) in (
      ("TRI", ctri, sorted_tri, TRI, 2, g2, ((:P,2,0),(:P,2,1),(:P⁻,1,2))),
      ("TET", ctet, sorted_tet, TET, 3, g3, ((:P⁻,1,0),)),
    )
  G, sq = unvarianced_inv(g), sqrt(det(g))
  for (F,r,k) in cases
    Uc, Vc = spaces(cmodel, fe(p,F,r,k))
    U , V  = spaces(model , fe(p,F,r,k))
    dΩc = Measure(Triangulation(cmodel), 2r+2)
    dΩ  = Measure(Triangulation(model) , 2r+2)

    mass(u,v)  = ∫( alt_inner(u,v) )dΩc
    stiff(u,v) = ∫( alt_inner(𝑑(u),𝑑(v)) )dΩc

    if k == 0
      # A 0-form is a function, of metric product `uv` and squared gradient
      # `g^{ij}∂ᵢu ∂ⱼv`: the mass and stiffness forms of Laplace–Beltrami.
      @test assemble_matrix(mass, Uc, Vc) ≈
        assemble_matrix((u,v) -> ∫( from_0form(u)*from_0form(v)*sq )dΩ, U, V)
      @test assemble_matrix(stiff, Uc, Vc) ≈ assemble_matrix(
        (u,v) -> ∫( (from_1form(𝑑(u))⋅(G⋅from_1form(𝑑(v))))*sq )dΩ, U, V)
    elseif k == 1
      @test assemble_matrix(mass, Uc, Vc) ≈
        assemble_matrix((u,v) -> ∫( (from_1form(u)⋅(G⋅from_1form(v)))*sq )dΩ, U, V)
    end

    if k == D
      # Raising the `D` indices of a top form divides by `det g`, which the
      # volume element of the scalar it becomes multiplies back by `√|det g|`.
      @test assemble_matrix(mass, Uc, Vc) ≈
        assemble_matrix((u,v) -> ∫( flat_vol_coeff(u)*flat_vol_coeff(v)/sq )dΩ, U, V)
    elseif k == D-1
      @test assemble_matrix(stiff, Uc, Vc) ≈ assemble_matrix(
        (u,v) -> ∫( flat_vol_coeff(𝑑(u))*flat_vol_coeff(𝑑(v))/sq )dΩ, U, V)
    end
  end
end

# ── The euclidean spelling is not the chart's ─────────────────────────────────

# `flat_vol_coeff(a ∧ ⋆ᴱb)` is `flat_alt_inner(a,b)`, a scalar: it takes the
# metric volume element of the measure all the same, and the metric the chart's
# `⋆` would have carried is the one missing. `LinearAlgebra.I` asks for that
# same euclidean operator.
let dΩ = Measure(Ωc,6), (U,V) = spaces(ctri, fe(TRI,:P⁻,1,1))
  euc = assemble_matrix(
    (u,v) -> ∫( flat_vol_coeff(𝑑(u) ∧ flat_hodge_star(𝑑(v))) )dΩ, U, V)
  @test euc ≈ assemble_matrix((u,v) -> ∫( alt_inner(𝑑(u),𝑑(v),LinearAlgebra.I) )dΩ, U, V)
  @test det(g2) * assemble_matrix((u,v) -> ∫( alt_inner(𝑑(u),𝑑(v)) )dΩ, U, V) ≈ euc
end

# ── The metric factor is read at the quadrature point ─────────────────────────

# `√|det g| = 1 + x¹` here, so a factor taken once per cell would not do.
vmodel = chart_model(sorted_tri, x -> SymTensorValue{Cov(2)}((1.0+x[1])^2,0.0,1.0))
@testset "a metric varying inside the cell, $F$r Λ$k" for (F,r,k) in
    ((:P,2,0),(:P⁻,1,1),(:P,2,1),(:P⁻,1,2))
  U, V = spaces(vmodel, fe(TRI,F,r,k))
  dΩv = Measure(Triangulation(vmodel), 2r+4)
  M = assemble_matrix((u,v) -> ∫( alt_inner(u,v) )dΩv, U, V)
  @test M ≈ assemble_matrix((u,v) -> ∫( u ∧ hodge_star(v) )dΩv, U, V)

  if k == 0
    Us, Vs = spaces(sorted_tri, fe(TRI,F,r,k))
    dΩ = Measure(Ω, 2r+4)
    G  = CellField(x -> TensorValue(1/(1.0+x[1])^2,0.0,0.0,1.0), Ω, PhysicalDomain())
    sq = CellField(x -> 1.0 + x[1], Ω, PhysicalDomain())
    @test assemble_matrix((u,v) -> ∫( alt_inner(𝑑(u),𝑑(v)) )dΩv, U, V) ≈ assemble_matrix(
      (u,v) -> ∫( (from_1form(𝑑(u))⋅(G⋅from_1form(𝑑(v))))*sq )dΩ, Us, Vs)
  end
end

# ── The exterior derivative takes no metric ───────────────────────────────────

# `𝑑` is the one operator of the calculus that does not read `g`, so an FE
# function of a chart model and its twin on the stock model share it.
let reffe = fe(TRI,:P,2,1)
  Vc, Vs = FESpace(ctri,reffe), FESpace(sorted_tri,reffe)
  dofs = collect(1.0:num_free_dofs(Vs))
  dc = 𝑑(FEFunction(Vc,dofs))(get_cell_points(Ωc))
  ds = 𝑑(FEFunction(Vs,dofs))(get_cell_points(Ω))
  @test all(cell -> dc[cell] ≈ ds[cell], eachindex(dc,ds))
end

# ── The metric mass matrix reproduces a form of the space ─────────────────────

# The `alt_inner` projection is an L² projection for the metric of the chart, so
# it is exact on the forms the space holds.
let dΩ = Measure(Ωc,6), (U,V) = spaces(ctri, fe(TRI,:P,2,1))
  f  = CellField(x -> ExteriorFormValue{1,Cov(2)}(1.0 + 2x[1] - x[2], 3.0 - x[1]),
                 Ωc, PhysicalDomain())
  ph = solve(AffineFEOperator((u,v) -> ∫( alt_inner(u,v) )dΩ,
                              v     -> ∫( alt_inner(f,v) )dΩ, U, V))
  @test sqrt(sum(∫( alt_inner(ph-f,ph-f) )dΩ)) < 1.0e-12
end

# ── A facet of a chart model carries no metric of its own ─────────────────────

# The induced metric is a different object, which the model does not hold, so
# assembling on a facet fails rather than read the cell metric.
let dΓ = Measure(BoundaryTriangulation(ctri),2), (U,V) = spaces(ctri, fe(TRI,:P,2,0))
  @test_throws ErrorException assemble_matrix((u,v) -> ∫( alt_inner(u,v) )dΓ, U, V)
end

end # module
