module MaterialHodgeStarsTests

using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Fields
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.TensorValues
using FillArrays
using Test
using LinearAlgebra

# A material (constitutive) Hodge star is the metric Hodge star composed with an
# endomorphism, on the `Λ^{D-K}` side (post) or on the `Λ^K` side (pre). The
# material is the plain matrix of that endomorphism in the stored components.

# The values a cell field takes at the first cell point of each of its cells.
function first_point_values(cf)
  a = cf(get_cell_points(get_triangulation(cf)))
  c = array_cache(a)
  [first(getindex!(c,a,cell)) for cell in 1:length(a)]
end

# The stored components of a form, and the form made of such a vector: what a
# material acts on, read and written without going through a contraction.
comps(ω) = collect(Tuple(ω))
form(::Val{J},::Val{D},c) where {J,D} = ExteriorFormValue{J,Cov(D)}(Tuple(c))

# The matrix of a linear operator on `Λ^J`, measured on the basis forms.
function operator_matrix(op,::Val{J},::Val{D}) where {J,D}
  L = binomial(D,J)
  hcat([comps(op(form(Val(J),Val(D),Float64.((1:L) .== i)))) for i in 1:L]...)
end

# ── The material is the matrix of the endomorphism, nothing more ──────────────

model3 = CartesianDiscreteModel((0.,1.,0.,1.,0.,1.),(2,2,2))
Ω3 = Triangulation(model3)
ω3 = CellField(x -> ExteriorFormValue{1,Cov(3)}(1.0 + x[1], 2.0, 3.0 + x[2]), Ω3, PhysicalDomain())
η3 = CellField(x -> ExteriorFormValue{2,Cov(3)}(1.0 + x[2], 2.0, 3.0), Ω3, PhysicalDomain())

# `K = 1` in dimension 3 has both sides act on `Λ²`, where the stored component
# contraction carries a `2!` that the material must not be charged with: the
# identity matrix gives `⋆` back, and a non-symmetric material acts as that
# matrix does on the stored components, in that order.
Nm = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 10.0]

@test all(first_point_values(postcomposed_material_hodge_star(Matrix(1.0I,3,3),1,Ω3)(ω3)) .≈
          first_point_values(hodge_star(ω3)))

⋆ᴺ = postcomposed_material_hodge_star(Nm,1,Ω3)
@test all(first_point_values(⋆ᴺ(ω3)) .≈
          map(s -> form(Val(2),Val(3),Nm*comps(s)), first_point_values(hodge_star(ω3))))

# the pre-composed material acts on `Λ^K`, so `K = 2` puts it on `Λ²` too
⋆ᴹ = precomposed_material_hodge_star(Nm,2,Ω3)
@test all(first_point_values(⋆ᴹ(η3)) .≈
          map(e -> flat_hodge_star(form(Val(2),Val(3),Nm*comps(e))), first_point_values(η3)))

# ── The two sides are conjugate: `N ∘ ⋆ = ⋆ ∘ M`, `M = ⋆⁻¹ ∘ N ∘ ⋆` ───────────

S3 = operator_matrix(flat_hodge_star,Val(1),Val(3))
@test all(first_point_values(precomposed_material_hodge_star(S3\(Nm*S3),1,Ω3)(ω3)) .≈
          first_point_values(⋆ᴺ(ω3)))

# ── Inverse: the side the material composes on swaps ──────────────────────────

# `inv(N ∘ ⋆) = (-1)^(K*(D-K)) ⋆ ∘ N⁻¹`, a pre-composed star of degree `D-K`.
# Dimension 2 is where that sign is not 1.
model2 = CartesianDiscreteModel((0.,1.,0.,1.),(2,2))
Ω2 = Triangulation(model2)
ω2 = CellField(x -> ExteriorFormValue{1,Cov(2)}(1.0 + 3x[1]*x[2], 2.0 + x[1]^2), Ω2, PhysicalDomain())
N2 = [2.0 1.0; 0.5 3.0]

⋆ᴺ₂ = postcomposed_material_hodge_star(N2,1,Ω2)
⋆ᴺ⁻¹ = precomposed_material_hodge_star(-inv(N2),1,Ω2)
@test all(first_point_values(⋆ᴺ⁻¹(⋆ᴺ₂(ω2))) .≈ first_point_values(ω2))

@test all(first_point_values(inv(⋆ᴺ₂)(⋆ᴺ₂(ω2))) .≈ first_point_values(ω2))
@test all(first_point_values(inv(⋆ᴺ₂)(⋆ᴺ₂(ω2))) .≈ first_point_values(⋆ᴺ⁻¹(⋆ᴺ₂(ω2))))

# `K = 1` in dimension 3 inverts into a pre-composed star of the 2-forms, whose
# material is the 3×3 `inv` of the original's
invᴺ₃ = inv(⋆ᴺ)
@test invᴺ₃ isa CellData.MaterialHodgeStar{2,3,:pre}
@test (⋆ᴺ).inverse == false && invᴺ₃.inverse == true
@test all(first_point_values(invᴺ₃(⋆ᴺ(ω3))) .≈ first_point_values(ω3))
@test all(first_point_values(⋆ᴺ(invᴺ₃(η3))) .≈ first_point_values(η3))

# and the pre-composed side inverts into the post-composed one
invᴹ₃ = inverse_map(⋆ᴹ)
@test invᴹ₃ isa CellData.MaterialHodgeStar{1,3,:post}
@test all(first_point_values(invᴹ₃(⋆ᴹ(η3))) .≈ first_point_values(η3))

# inverting twice gives the operator back, the material field included
@test inv(invᴺ₃) isa typeof(⋆ᴺ)
@test inv(invᴺ₃).inverse == (⋆ᴺ).inverse
@test inv(invᴺ₃).m === (⋆ᴺ).m

# ── The kinds a material may be given as ──────────────────────────────────────

# in dimension 3, where the `2!` of `Λ²` tells a material that is rescaled once
# at construction from one that is rescaled at each point, or not at all
N3v = TensorValue{3,3}(Tuple(Nm))
ref3 = first_point_values(⋆ᴺ(ω3))
for material in (N3v, Nm, Fill(N3v,num_cells(model3)), x -> N3v,
                 CellField(x -> N3v, Ω3, PhysicalDomain()))
  @test all(first_point_values(postcomposed_material_hodge_star(material,1,Ω3)(ω3)) .≈ ref3)
end

# a material of the wrong size, a degree out of range, and a material that is
# not the matrix of an endomorphism
@test_throws AssertionError postcomposed_material_hodge_star(N2,1,Ω3)
@test_throws AssertionError postcomposed_material_hodge_star(x -> N3v,4,Ω3)
@test_throws AssertionError postcomposed_material_hodge_star(x -> VectorValue(1.0,2.0),1,Ω2)

# ── An FE basis goes through either side alike ────────────────────────────────

V = FESpace(model2,ReferenceFE(QUAD,:Q⁻,1,1,Float64),conformity=:HΛ)
dΩ2 = Measure(Ω2,2)
S2 = operator_matrix(flat_hodge_star,Val(1),Val(2))
⋆ᴹ₂ = precomposed_material_hodge_star(S2\(N2*S2),1,Ω2)

@test assemble_matrix((u,v) -> ∫( v ∧ ⋆ᴺ₂(u) )dΩ2, V, V) ≈
      assemble_matrix((u,v) -> ∫( v ∧ ⋆ᴹ₂(u) )dΩ2, V, V)

# ── Against the wedge, on a chart of non-trivial metric ───────────────────────

const g = SymTensorValue{Cov(2)}(2.0,0.5,1.0)

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

cmodel = ChartModel(model2,Fill(GenericField(x -> g),num_cells(model2)))
Ωc  = Triangulation(cmodel)
dΩc = Measure(Ωc,2)
ωc  = CellField(x -> ExteriorFormValue{1,Cov(2)}(1.0 + 3x[1]*x[2], 2.0 + x[1]^2), Ωc, PhysicalDomain())
ηc  = CellField(x -> ExteriorFormValue{1,Cov(2)}(0.5, 1.0), Ωc, PhysicalDomain())

# The material weights the pointwise pairing of the forms and leaves the volume
# form alone. The left hand sides are top forms, which integrate with no metric
# weight, so the `√|det g|` of `μ` is the one the scalars on the right are
# integrated against.
md  = compute_metric_data(g)
⋆ᴺᶜ = postcomposed_material_hodge_star(N2,1,Ωc)

# `ω ∧ ⋆ᴺ(η) = (⋆ω ⨟ N(⋆η)) μ`: `N` weights the pairing of the `D-K`-forms,
# which is the space it is an endomorphism of, with no conjugation anywhere.
sη   = hodge_star(ExteriorFormValue{1,Cov(2)}(0.5, 1.0), md)
Nsηc = CellField(x -> form(Val(1),Val(2),N2*comps(sη)), Ωc, PhysicalDomain())
@test sum(∫( ωc ∧ ⋆ᴺᶜ(ηc) )dΩc) ≈ sum(∫( hodge_star(ωc) ⨟ Nsηc )dΩc)

# `ω ∧ ⋆ᴹ(η) = (ω ⨟ M(η)) μ` is the same identity read on the `K`-forms, `⋆`
# being an isometry of `⨟`; equating the two gives `M = ⋆⁻¹ ∘ N ∘ ⋆` back.
Sc  = operator_matrix(w -> hodge_star(w,md),Val(1),Val(2))
Mηc = CellField(x -> form(Val(1),Val(2),(Sc\(N2*Sc))*[0.5,1.0]), Ωc, PhysicalDomain())
@test sum(∫( ωc ∧ ⋆ᴺᶜ(ηc) )dΩc) ≈ sum(∫( ωc ⨟ Mηc )dΩc)

# the `⋆` of the chart metric is the one composed, not the Euclidean one
η2 = CellField(x -> ExteriorFormValue{1,Cov(2)}(0.5, 1.0), Ω2, PhysicalDomain())
@test !(sum(∫( ωc ∧ ⋆ᴺᶜ(ηc) )dΩc) ≈ sum(∫( ω2 ∧ ⋆ᴺ₂(η2) )dΩ2))

@test all(first_point_values(inv(⋆ᴺᶜ)(⋆ᴺᶜ(ηc))) .≈ first_point_values(ηc))

# ── Differentiation: the operator commutes with ∂ as the ⋆ it composes does ───

# The material contracts the form linearly, so a material Hodge star is a
# bilinear operation of its material and its argument, and the Leibniz rule
# commutes it with `∂` at the orders the composed Hodge star is derived at. The
# reference is the same map written in closed form, which needs no such rule.
ω2fun(x) = ExteriorFormValue{1,Cov(2)}(1.0 + 3x[1]*x[2]^2, 2.0 + x[1]^2)
ω2d = CellField(ω2fun, Ω2, PhysicalDomain())
refᴺ = CellField(x -> form(Val(1),Val(2),N2*comps(flat_hodge_star(ω2fun(x)))), Ω2, PhysicalDomain())

# the pre-composed conjugate is the same map, so it has the same derivatives
for material_star in (⋆ᴺ₂, ⋆ᴹ₂),
    k in (exterior_derivative, partial_derivatives, second_partial_derivatives)
  @test all(first_point_values(k(material_star(ω2d))) .≈ first_point_values(k(refᴺ)))
end

# On a chart, the metric Hodge star is derived at order 1 only, and so is a
# material star composing it.
ωcd = CellField(ω2fun, Ωc, PhysicalDomain())
xc  = get_cell_points(Ωc)
refᴺᶜ = CellField(x -> form(Val(1),Val(2),N2*comps(hodge_star(ω2fun(x),md))), Ωc, PhysicalDomain())
for k in (exterior_derivative, partial_derivatives)
  @test all(collect(k(⋆ᴺᶜ(ωcd))(xc)) .≈ collect(k(refᴺᶜ)(xc)))
end
@test_throws ErrorException collect(second_partial_derivatives(⋆ᴺᶜ(ωcd))(xc))
@test_throws ErrorException collect(second_partial_derivatives(hodge_star(ωcd))(xc))

# ── An FE basis is derived through either side alike ──────────────────────────

# `𝑑(⋆ᴺu)` of an FE basis goes through the lazy path, where the material star
# carries the rule of a bilinear operation. The material mixes the components of
# `⋆u`, so this does not vanish on a space of divergence free Whitney forms, as
# `𝑑(⋆u)` itself does.
aᴺ(u,v) = ∫( exterior_derivative(⋆ᴺ₂(u)) ⨟ exterior_derivative(⋆ᴺ₂(v)) )dΩ2
aᴹ(u,v) = ∫( exterior_derivative(⋆ᴹ₂(u)) ⨟ exterior_derivative(⋆ᴹ₂(v)) )dΩ2
Aᴺ = assemble_matrix(aᴺ, V, V)
@test Aᴺ ≈ assemble_matrix(aᴹ, V, V)
@test norm(Aᴺ) > 0
@test iszero(assemble_matrix((u,v) -> ∫(
  exterior_derivative(hodge_star(u)) ⨟ exterior_derivative(hodge_star(v)) )dΩ2, V, V))

end # module
