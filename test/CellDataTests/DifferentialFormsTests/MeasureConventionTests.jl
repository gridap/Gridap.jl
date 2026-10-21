module MeasureConventionTests

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

# The volume element of `∫( … )dΩ`: `meas(Jt)·√|det g|` of the model's metric on
# a density, and the euclidean element on a differential form, which carries the
# volume factor of its own Hodge star.

# ── Chart models of given metric ──────────────────────────────────────────────

struct CellMetric{F} <: Field
  g::F
end
Gridap.Arrays.evaluate!(cache, m::CellMetric, x::Point) = m.g(x)

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

chart_model(model,g) = ChartModel(model,Fill(CellMetric(g),num_cells(model)))

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(4,4))
Ω  = Triangulation(model)
dΩ = Measure(Ω,3)

# g = [[2, 1/2], [1/2, 1]], covariant as a metric must be
const g = SymTensorValue{Cov(2)}(2.0,0.5,1.0)
cmodel = chart_model(model, x -> g)
Ωc  = Triangulation(cmodel)
dΩc = Measure(Ωc,3)
sqrt_g = sqrt(det(g))

# ── A density is integrated against the metric volume element ─────────────────

@test sum(∫( 1.0 )dΩ)  ≈ 1.0
@test sum(∫( 1.0 )dΩc) ≈ sqrt_g

# The factor is that of the point it multiplies, not of the cell: `√|det g|` is
# `1 + x¹` here, of integral 3/2 over the unit square.
vmodel = chart_model(model, x -> SymTensorValue{Cov(2)}((1.0 + x[1])^2,0.0,1.0))
dΩv = Measure(Triangulation(vmodel),3)
@test sum(∫( 1.0 )dΩv) ≈ 3/2

# ── A differential form is integrated as the coefficient of a top form ────────

ω2(x) = ExteriorFormValue{2,Cov(2)}(1.0 + x[1]*x[2])
ω1(x) = ExteriorFormValue{1,Cov(2)}(1.0 + 3x[1]*x[2], 2.0 + x[1]^2)
ω0(x) = ExteriorFormValue{0,2}(x[1])

ω2_cf , ω1_cf , ω0_cf  = map(f -> CellField(f,Ω ,PhysicalDomain()), (ω2,ω1,ω0))
ω2c_cf, ω1c_cf, ω0c_cf = map(f -> CellField(f,Ωc,PhysicalDomain()), (ω2,ω1,ω0))

# A top form is its single coefficient, and takes no metric volume factor
@test sum(∫( ω2_cf  )dΩ)  ≈ sum(∫( flat_vol_coeff(ω2_cf) )dΩ)
@test sum(∫( ω2c_cf )dΩc) ≈ sum(∫( flat_vol_coeff(ω2_cf) )dΩ)

# Reducing it by hand gives a density instead, which does take the factor: the
# value type of the integrand is what the convention reads.
@test sum(∫( flat_vol_coeff(ω2c_cf) )dΩc) ≈ sqrt_g * sum(∫( flat_vol_coeff(ω2_cf) )dΩ)

# A 0-form is a function, and takes the factor
@test sum(∫( ω0_cf  )dΩ)  ≈ 1/2
@test sum(∫( ω0c_cf )dΩc) ≈ sqrt_g / 2

# The components of a form of any other degree depend on the basis they are
# written in, so it has no integral.
@test_throws ErrorException sum(∫( ω1_cf  )dΩ)
@test_throws ErrorException sum(∫( ω1c_cf )dΩc)

# A top form of the model is not a top form of its boundary, which needs the
# form pulled back to it; a 0-form is a function on either.
Γ  = BoundaryTriangulation(model)
dΓ = Measure(Γ,3)
@test_throws ErrorException sum(∫( CellField(ω2,Γ,PhysicalDomain()) )dΓ)
@test sum(∫( CellField(ω0,Γ,PhysicalDomain()) )dΓ) ≈ 2.0

# The volume element of a facet of a chart model is that of the metric induced
# on it, which the model does not carry, so integrating on one fails rather than
# take the euclidean element.
@test_throws ErrorException sum(∫( 1.0 )Measure(BoundaryTriangulation(cmodel),3))

# ── `∫(a ∧ ⋆b)dΩ == ∫(alt_inner(a,b))dΩ` for every degree and dimension ───────

# The convention is one rule read from both ends, `a ∧ ⋆b = √|det g| (a|b) dx¹∧…∧dx^D`:
# the metric enters the form spelling through `⋆` and the scalar one through the
# volume element. The tensor inner product differs from both by `K!`.
const g3 = SymTensorValue{Cov(3)}(2.0,0.5,0.0,1.0,0.25,1.5)
model3  = CartesianDiscreteModel((0.0,1.0,0.0,1.0,0.0,1.0),(2,2,2))
cmodel3 = chart_model(model3, x -> g3)
Ωc3  = Triangulation(cmodel3)
dΩc3 = Measure(Ωc3,3)

for (D,Ωk,dΩk) in ((2,Ωc,dΩc), (3,Ωc3,dΩc3))
  for K in 0:D
    L = binomial(D,K)
    af(x) = ExteriorFormValue{K,Cov(D)}(ntuple(i -> 1.0 + i*x[1], L)...)
    bf(x) = ExteriorFormValue{K,Cov(D)}(ntuple(i -> 2.0 - i*x[2], L)...)
    a, b = map(f -> CellField(f,Ωk,PhysicalDomain()), (af,bf))
    @test sum(∫( a ∧ hodge_star(b) )dΩk) ≈ sum(∫( alt_inner(a,b) )dΩk)
    @test sum(∫( inner(a,b) )dΩk) ≈ factorial(K) * sum(∫( alt_inner(a,b) )dΩk)
    # The `⋆` of a top form divides its coefficient by `√|det g|`, which the
    # volume element of the 0-form it becomes restores.
    @test sum(∫( hodge_star(a ∧ hodge_star(b)) )dΩk) ≈ sum(∫( alt_inner(a,b) )dΩk)
  end
end

# ── Assembly on a form valued space takes no `flat_vol_coeff` by hand ─────────

V0 = FESpace(model,ReferenceFE(QUAD,lagrangian,ExteriorFormValue{0,2,Float64},1))
U0 = TrialFESpace(V0)

A    = assemble_matrix((u,v) -> ∫( 𝑑(u) ∧ hodge_star(𝑑(v)) )dΩ, U0, V0)
A_vc = assemble_matrix((u,v) -> ∫( flat_vol_coeff(𝑑(u) ∧ flat_hodge_star(𝑑(v))) )dΩ, U0, V0)
@test A ≈ A_vc

f = CellField(x -> ExteriorFormValue{0,2}(2π^2*sin(π*x[1])*sin(π*x[2])),Ω,PhysicalDomain())
b    = assemble_vector(v -> ∫( f * v )dΩ, V0)
b_hs = assemble_vector(v -> ∫( f ∧ hodge_star(v) )dΩ, V0)
@test b ≈ b_hs

# The stiffness matrix of the 0-form space is that of the Lagrangian space of
# the same nodes, whose dofs it shares.
Vs = FESpace(model,ReferenceFE(QUAD,lagrangian,Float64,1))
Us = TrialFESpace(Vs)
@test A ≈ assemble_matrix((u,v) -> ∫( ∇(u)⋅∇(v) )dΩ, Us, Vs)

# ── A reference domain integral takes no metric factor ────────────────────────

# The metric factor accompanies `meas(Jt)`: it is the chart's half of the change
# of variables, and an integral that does not change variables takes neither.
quad_ref = CellQuadrature(Ωc,3;integration_domain_style=ReferenceDomain())
@test sum(∫( 1.0 )quad_ref) ≈ num_cells(cmodel)

# A Dirac delta integrates in the physical domain, weights included, and is a
# measure of its own: it takes no volume element at all.
δ = DiracDelta(cmodel,Point(0.5,0.25))
@test sum(δ(CellField(x -> x[1],Ωc,PhysicalDomain()))) ≈ 0.5

# ── FE assembly on a stock model is untouched ─────────────────────────────────

# The Galerkin solution of `−Δu = 1` with `u = 0` on the boundary satisfies
# `a(uh,uh) = l(uh)`, which holds of the assembled system alone.
V = FESpace(model,ReferenceFE(lagrangian,Float64,1),dirichlet_tags="boundary")
U = TrialFESpace(V,0.0)
uh = solve(AffineFEOperator((u,v) -> ∫( ∇(u)⋅∇(v) )dΩ, v -> ∫( v )dΩ, U, V))
@test sum(∫( ∇(uh)⋅∇(uh) )dΩ) ≈ sum(∫( uh )dΩ)
@test sum(∫( uh )dΩ) > 0

end # module
