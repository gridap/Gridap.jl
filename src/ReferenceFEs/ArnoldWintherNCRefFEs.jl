"""
    struct ArnoldWintherNC <: ReferenceFEName end

Reference FE name for the nonconforming Arnold--Winther stress element. See
[`ArnoldWintherNCRefFE`](@ref). Its singleton instance is [`aw_nc`](@ref).
"""
struct ArnoldWintherNC <: ReferenceFEName end

"""
    const aw_nc = ArnoldWintherNC()

Singleton of the [`ArnoldWintherNC`](@ref) reference FE name.
"""
const aw_nc = ArnoldWintherNC()

Pushforward(::Type{ArnoldWintherNC}) = DoubleContraVariantPiolaMap()

"""
    ArnoldWintherNCRefFE(::Type{T}, p::Polytope{2})

The nonconforming Arnold--Winther reference FE on the triangle `K`, with `T` the
scalar type: 15 DoFs and, writing `S` for the symmetric 2×2 matrices,

    AWnc(K) = {τ ∈ P₂(K;S) : (n⋅τn)|ₑ ∈ P₁(e) ∀ e ⊂ ∂K}

of dimension 15 [Arnold & Winther, M3AS 13 (2003) 295].

Assembled as the augmented element of [Kirby, SMAI-JCM 4 (2018) 197, §5.5]: the
DoFs *and* the constraint functionals are declared together over an ambient
prebasis they form a dual basis for, the generalized Vandermonde is inverted, and
the columns dual to the DoFs are kept — they lie in the constrained space, being
annihilated by every constraint. Note `length(get_prebasis(reffe))` is therefore
larger than `num_dofs(reffe)`.

## Prebasis

`P₂(K;S)`, of dimension 18 — the ambient space, not `AWnc(K)` itself.

## Moments

Per edge `e`, with unit tangent `t`, normal `n = R t` and `μᵢ` the
L²(e)-orthonormal Legendre basis:

    ℓ^{nn,i}_e(τ) = ∫ₑ (n⋅τn) μᵢ ds,    i = 0, 1
    ℓ^{nt,i}_e(τ) = ∫ₑ (n⋅τt) μᵢ ds,    i = 0, 1

four per edge; and over the cell the three independent components,

    ℓ^{K,c}(τ) = ∫_K τ ⊙ E_c dK,    E_c ∈ {e₁⊗e₁, e₁⊗e₂, e₂⊗e₂}

`4·3 + 3 = 15`.

## Constraints

    ∫ₑ (n⋅τn) μ ds = 0    ∀ μ ∈ P₂(e) ∩ P₁(e)^⊥,  ∀ e       (1 per edge)

Three in all: `(n⋅τn)|ₑ` lies in `P₂(e)` a priori and is pinned to its `P₁` part.
`15 + 3 = 18`.
"""
function ArnoldWintherNCRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The nonconforming Arnold-Winther element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), SymTensorValue{2,T}, 2)
  fb = LegendreBasis(Val(1), T, 1)                            # μ₀, μ₁ : the DoF weights
  # μ₂ : the constraint, the degree-2 Legendre polynomial alone. It spans
  # P₂(e) ∩ P₁(e)^⊥, so the single moment against it states (n⋅τn)|ₑ ∈ P₁(e).
  gb = LegendreBasis(Val(1), T, 2, _pk_minus_pq_filter(2, 1))
  cb = MonomialBasis(Val(2), T, 0, _p_filter)                 # the constant on the cell
  Ei = (
    ConstantField(TensorValue(one(T), zero(T), zero(T), zero(T))),   # e₁⊗e₁
    ConstantField(TensorValue(zero(T), zero(T), one(T), zero(T))),   # e₁⊗e₂
    ConstantField(TensorValue(zero(T), zero(T), zero(T), one(T))),   # e₂⊗e₂
  )

  function nnmom(φ, μ, ds)
    n = _edge_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ, n)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φn, n), μ)
  end
  function ntmom(φ, μ, ds)
    n = _edge_normal(ds)
    t = get_edge_tangent(ds)
    φn = Broadcasting(Operation(⋅))(φ, n)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φn, t), μ)
  end
  cmom(E) = (φ, μ, ds) -> Broadcasting(Operation(*))(   # ∫_K τ ⊙ E dK
    Broadcasting(Operation(⊙))(φ, E), μ
  )

  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  moments = Tuple[
    (edges, nnmom, fb), (edges, ntmom, fb), (edges, nnmom, gb), # Edge moments
    [(cell, cmom(E), cb) for E in Ei]...,                       # Cell moments
  ]

  # Note the constraint weight is degree 2 while the DoF weights are degree 1, so
  # the three moment entries on an edge use three different quadratures -- the
  # configuration that needs `MomentBasedDofBasis` to offset each entry's moment
  # rows by the face's already-used node count.
  full = MomentBasedDofBasis(p, prebasis, moments)
  nedges = num_faces(p, 1)
  dof_ids = vcat(reduce(vcat, [5*(k - 1) .+ (1:4) for k in 1:nedges]), 5*nedges .+ (1:3))
  ndofs = length(dof_ids)

  @check length(full) == length(prebasis) """\n
  The augmented element needs as many functionals as prebasis functions, got
  $(length(full)) and $(length(prebasis)).
  """
  shapefuns = linear_combination(inv(evaluate(full, prebasis))[:, dof_ids], prebasis)
  dofs = restrict(full, dof_ids)

  @check maximum(abs, evaluate(dofs, shapefuns) - Matrix{Float64}(I, ndofs, ndofs)) < 1e-10 """\n
  The augmented AWnc Vandermonde is singular or badly conditioned.
  """
  GenericRefFE{ArnoldWintherNC}(
    ndofs, p, prebasis, dofs, DivConformity(), nothing,
    get_face_own_moments(dofs), shapefuns
  )
end

function ReferenceFE(p::Polytope, ::ArnoldWintherNC, ::Type{T}) where T
  ArnoldWintherNCRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::ArnoldWintherNC, ::Type{T}, order) where T
  @notimplementedif order != 2 """\n
  The nonconforming Arnold-Winther element exists for order 2 only, got $order.
  """
  ArnoldWintherNCRefFE(T, p)
end

# Identity for every admissible vertex permutation of every face: the DoFs of an
# edge are ordered by moment kind and then by Legendre degree, which both
# adjacent cells agree on, and reversing the edge only changes the signs of the
# odd-degree ones, which the change of basis carries.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{ArnoldWintherNC}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end
