"""
    struct ArnoldWintherC <: ReferenceFEName end

Reference FE name for the conforming Arnold--Winther stress element. See
[`ArnoldWintherCRefFE`](@ref). Its singleton instance is [`aw_c`](@ref).
"""
struct ArnoldWintherC <: ReferenceFEName end

"""
    const aw_c = ArnoldWintherC()

Singleton of the [`ArnoldWintherC`](@ref) reference FE name.
"""
const aw_c = ArnoldWintherC()

Pushforward(::Type{ArnoldWintherC}) = DoubleContraVariantPiolaMap()

"""
    ArnoldWintherCRefFE(::Type{T}, p::Polytope{2})

The conforming Arnold--Winther reference FE on the triangle `K`, with `T` the
scalar type: 24 DoFs and, writing `S` for the symmetric 2×2 matrices,

    AWc(K) = {τ ∈ P₃(K;S) : div τ ∈ P₁(K;R²)}

of dimension 24 [Arnold & Winther, Numer. Math. 92 (2002) 401].

Assembled as the augmented element of [Kirby, SMAI-JCM 4 (2018) 197, §5.5], as
[`ArnoldWintherNCRefFE`](@ref) is — see there for what the construction does.

## Prebasis

`P₃(K;S)`, of dimension 30 — the ambient space, not `AWc(K)` itself.

## Moments

At each vertex `v` the three independent components,

    ℓ^{v,c}(τ) = τ(v) ⊙ E_c,    E_c ∈ {e₁⊗e₁, e₁⊗e₂, e₂⊗e₂}

per edge `e`, with unit tangent `t`, normal `n = R t` and `μᵢ` the
L²(e)-orthonormal Legendre basis,

    ℓ^{nn,i}_e(τ) = ∫ₑ (n⋅τn) μᵢ ds,    ℓ^{nt,i}_e(τ) = ∫ₑ (n⋅τt) μᵢ ds,   i = 0, 1

and over the cell the same three components,

    ℓ^{K,c}(τ) = ∫_K τ ⊙ E_c dK.

`3·3 + 4·3 + 3 = 24`.

The `E_c` are full matrices, not the component basis of `S`: contracting two
symmetric matrices sums over both off-diagonal slots, so the symmetric
off-diagonal basis element would give `2τ₁₂`. The vertex DoFs must be exactly
`(τ₁₁, τ₁₂, τ₂₂)`, since the vertex block of the change of basis is the matrix
of `H ↦ J H Jᵀ` written in those components.

## Constraints

    ∫_K (div τ)_c q dK = 0    ∀ q ∈ P₂(K) ∩ P₁(K)^⊥,  c = 1, 2      (3 each)

Six in all: each component of `div τ` lies in `P₂(K)` a priori and is pinned to
its `P₁` part. `24 + 6 = 30`.
"""
function ArnoldWintherCRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The conforming Arnold-Winther element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), SymTensorValue{2,T}, 3)

  # DoF moments
  vb = MonomialBasis(Val(0), T, 0)                # the constant at a vertex
  fb = LegendreBasis(Val(1), T, 1)                # μ₀, μ₁ on the edge
  cb = MonomialBasis(Val(2), T, 0, _p_filter)     # the constant on the cell
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
  Emom(E) = (φ, μ, ds) -> Broadcasting(Operation(*))(
    Broadcasting(Operation(⊙))(φ, E), μ
  )

  verts = get_dimrange(p, 0)
  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  moments = Tuple[
    [(verts, Emom(E), vb) for E in Ei]...,  # Vertex moments
    (edges, nnmom, fb), (edges, ntmom, fb), # Edge moments
    [(cell, Emom(E), cb) for E in Ei]...,   # Cell moments
  ]

  # Constraint moments.
  # P₂(K) ∩ P₁(K)^⊥: the degree-2 Dubiner polynomials, i.e. the members of an
  # L²(K)-orthonormal basis of P₂ that are orthogonal to every linear. A filtered
  # Legendre basis would not do -- it is not orthogonal on a simplex.
  qb = DubinerBasis(Val(2), T, 2, _p_complement_filter(1))
  Ej = (
    ConstantField(VectorValue(one(T), zero(T))),
    ConstantField(VectorValue(zero(T), one(T)))
  )
  divmom(e) = (φ, μ, ds) -> Broadcasting(Operation(*))(
    Broadcasting(Operation(⋅))(Broadcasting(Operation(tr))(φ), e), μ
  )
  constraints = [ (cell, divmom(e), qb) for e in Ej ]

  # the 24 DoFs all take the identity operator, the 6 divergence constraints take
  # ∇, so they are two bases joined with `vcat`.
  dofs = MomentBasedDofBasis(p, prebasis, moments)
  cons = MomentBasedDofBasis(p, prebasis, constraints, ∇)
  full = vcat(dofs, cons)
  ndofs = length(dofs)

  @check length(full) == length(prebasis) """\n
  The augmented element needs as many functionals as prebasis functions, got
  $(length(full)) and $(length(prebasis)).
  """
  # the DoFs come first, so no `restrict` is needed: the slice suffices
  shapefuns = linear_combination(inv(evaluate(full, prebasis))[:, 1:ndofs], prebasis)

  @check maximum(abs, evaluate(dofs, shapefuns) - Matrix{Float64}(I, ndofs, ndofs)) < 1e-10 """\n
  The augmented AWc Vandermonde is singular or badly conditioned.
  """
  GenericRefFE{ArnoldWintherC}(
    ndofs, p, prebasis, dofs, DivConformity(), nothing,
    get_face_own_moments(dofs), shapefuns
  )
end

function ReferenceFE(p::Polytope, ::ArnoldWintherC, ::Type{T}) where T
  ArnoldWintherCRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::ArnoldWintherC, ::Type{T}, order) where T
  @notimplementedif order != 3 """\n
  The conforming Arnold-Winther element exists for order 3 only, got $order.
  """
  ArnoldWintherCRefFE(T, p)
end

# Identity for every admissible vertex permutation of every face. The vertex DoFs
# are components in the global Cartesian frame, hence the same functionals for
# every cell touching the vertex; the edge DoFs are ordered by moment kind and
# Legendre degree, and reversing an edge only changes the signs of the odd-degree
# ones, which the change of basis carries.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{ArnoldWintherC}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end
