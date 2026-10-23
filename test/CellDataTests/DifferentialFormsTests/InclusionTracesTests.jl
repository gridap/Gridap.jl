module InclusionTracesTests

using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.TensorValues
using Test
using LinearAlgebra

using Gridap.CellData: induced_orientation

# The trace of a form is the pullback of the inclusion, `trΓ(ω) = ι*ω`. What
# fixes it is Stokes' theorem,
#
#   ∫_Ω 𝑑ω = ∮_{∂Ω} trΓ(ω),
#
# which is what the orientation and the measure of `Γ` are there for. Every
# integral of a top form below is written in both ways, as the integral of the
# form and as the integral of the function it stars to,
#
#   ∫( ω )dΓ = ∫( ⋆ω )dΓ,
#
# the two agreeing only if the `⋆` of `Γ` and its measure carry the same
# orientation.

# The values a cell field takes at the first cell point of each of its cells.
function first_point_values(cf)
  a = cf(get_cell_points(get_triangulation(cf)))
  c = array_cache(a)
  [first(getindex!(c,a,cell)) for cell in 1:length(a)]
end

test_form(::Val{2}) = x -> ExteriorFormValue{1,Cov(2)}(x[1]^2*x[2], 3.0*x[2]^2 + x[1])
test_form(::Val{3}) = x -> ExteriorFormValue{2,Cov(3)}(
  x[3]^2 + x[1]*x[2], x[1]*x[2], x[1]^2*x[2] + x[2]*x[3])

function stokes_terms(model, ::Val{D}, deg=6) where D
  Ω = Triangulation(model); dΩ = Measure(Ω,deg)
  Γ = BoundaryTriangulation(model); dΓ = Measure(Γ,deg)
  trΓ = get_trace_operator(Γ,Ω)
  ω = CellField(test_form(Val(D)), Ω, PhysicalDomain())
  (sum(∫( 𝑑(ω) )dΩ), sum(∫( trΓ(ω) )dΓ), sum(∫( ⋆(trΓ(ω)) )dΓ))
end

# ── Stokes' theorem, in both writings ─────────────────────────────────────────

# The simplex meshes matter here: `simplexify` emits cells of both
# orientations, which only a trace that carries the induced one survives.
for (D, dom, parts) in ((2, (0.,1.,0.,1.),          (2,2)),
                        (2, (0.,1.,0.,1.),          (3,3)),
                        (3, (0.,1.,0.,1.,0.,1.),    (1,1,1)),
                        (3, (0.,1.,0.,1.,0.,1.),    (2,2,2)))
  cartesian = CartesianDiscreteModel(dom,parts)
  for model in (cartesian, simplexify(cartesian))
    exterior_derivative, as_form, as_function = stokes_terms(model, Val(D))
    @test as_form     ≈ exterior_derivative
    @test as_function ≈ exterior_derivative
  end
end

# ── The trace keeps the tangential part and drops the rest ────────────────────

n = 4
model = CartesianDiscreteModel((0.,1.,0.,1.),(n,n))
Ω = Triangulation(model)

# tag 7 is an edge of constant `x¹`, so `dx¹` pulls back to zero there and
# `dx²` to the one component of a top form of that edge. That component is in
# the reference coordinates of the edge, where the tangent has length `h`.
Γ7  = BoundaryTriangulation(model, tags=[7])
tr7 = get_trace_operator(Γ7,Ω)
dx1 = CellField(ExteriorFormValue{1,Cov(2)}(1.0,0.0), Ω, PhysicalDomain())
dx2 = CellField(ExteriorFormValue{1,Cov(2)}(0.0,1.0), Ω, PhysicalDomain())

@test all(v -> Tuple(v) == (0.0,), first_point_values(tr7(dx1)))
@test all(v -> Tuple(v) == (1/n,), first_point_values(tr7(dx2)))

# a trace needs a triangulation of lower dimension to be taken on
@test_throws AssertionError get_trace_operator(Ω,Ω)

# ── The metric and the orientation Γ inherits ─────────────────────────────────

Γ = BoundaryTriangulation(model)
g, g_inv, vol = get_metric_data(Γ)

# the edges are parametrised by the unit segment, so they inherit `h²`
@test all(v -> Tuple(v) == (1/n^2,), first_point_values(g))
@test all(v -> Tuple(v) == (n^2*1.0,), first_point_values(g_inv))

# the volume element is their length, and carries no orientation: a metric has
# none. The sign that makes the two writings of an integral agree is held by the
# triangulation, and reaches the measure and the `⋆` from there.
@test all(v -> v ≈ 1/n, first_point_values(vol))

ε = induced_orientation(Γ)
@test length(ε) == num_cells(Γ)
@test sort(unique(ε)) == [-1.0, 1.0]

# ── A Neumann term assembles the same in every writing ────────────────────────

V   = FESpace(model, ReferenceFE(QUAD,:Q⁻,1,1,Float64), conformity=:HΛ)
Γ8  = BoundaryTriangulation(model, tags=[8])
dΓ8 = Measure(Γ8,4)
tr8 = get_trace_operator(Γ8,Ω)
p̄   = CellField(ExteriorFormValue{0,Cov(1)}(-1.0), Γ8, PhysicalDomain())

b_form = assemble_vector(v -> ∫( p̄ ∧ tr8(v)    )dΓ8, V)
b_star = assemble_vector(v -> ∫( p̄ ∧ ⋆(tr8(v)) )dΓ8, V)

@test b_form ≈ b_star
@test !iszero(b_form)

# ── An interior facet orients its `⋆` the same way ────────────────────────────

# `jump` and `mean` are taken on `Ω` and the trace after them, a lifted operator
# bringing its arguments to a common domain as it is built. `Λ` is oriented by
# its `plus` side, so the two writings agree there as they do on a boundary.
Λ   = SkeletonTriangulation(model)
dΛ  = Measure(Λ,4)
trΛ = get_trace_operator(Λ,Ω)
ωΛ  = CellField(test_form(Val(2)), Ω, PhysicalDomain())

@test sum(∫( trΛ(mean(ωΛ)) )dΛ) ≈ sum(∫( ⋆(trΛ(mean(ωΛ))) )dΛ)
@test !iszero(sum(∫( trΛ(mean(ωΛ)) )dΛ))

end # module
