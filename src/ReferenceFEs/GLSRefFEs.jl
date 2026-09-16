"""
    struct GopalakrishnanLedererSchoberl <: ReferenceFEName end

Reference FE name for the Gopalakrishnan--Lederer--Schoberl element (second
kind), a normal-tangential-continuous traceless-matrix-valued element for the
mass-conserving mixed stress formulation of Stokes flow. See [`GLSRefFE`](@ref).
Its singleton instance is [`gls`](@ref).
"""
struct GopalakrishnanLedererSchoberl <: ReferenceFEName end

"""
    const gls = GopalakrishnanLedererSchoberl()

Singleton of the [`GopalakrishnanLedererSchoberl`](@ref) reference FE name.
"""
const gls = GopalakrishnanLedererSchoberl()

Pushforward(::Type{GopalakrishnanLedererSchoberl}) = CoContraVariantPiolaMap()

"""
    GLSRefFE(::Type{T}, K::Polytope{2}, order::Integer)

The Gopalakrishnan--Lederer--Schoberl reference FE of the second kind, degree
`r = order` on the triangle `K`, with `T` the scalar type. Writing `M₀` for the
traceless 2×2 tensors,

    GLS_r(K) = P_r(K;M₀)

of dimension `3(r+1)(r+2)/2`, for any `r ≥ 0`
[Gopalakrishnan, Lederer & Schoberl, SIAM J. Numer. Anal. 58 (2020) 706].

This element is divergence conforming in the sense that the normal-normal trace
on facets `n⋅τ⋅n` is pointwise continuous.

The implementation follows the augmented element approach of [Kirby, SMAI-JCM 4 (2018) 197].

# Extended help

## Prebasis

The prebasis is taken as `P_r(K;M)`, of dimension `4(r+1)(r+2)/2`, with constraints

    ∫_K tr(M) q dK = 0    ∀ q ∈ P_r(K)

enforced using moments, yielding the `3(r+1)(r+2)/2` DoFs of the GLS element.

## Moments

Per edge `e` with unit tangent `t`, normal `n = R t` and `μᵢ` the
L²(e)-orthonormal Legendre basis of `P_r(e)`; and over the cell, the same
normal-tangential functional of each edge's frame against a basis `{q}` of
`P_{r-1}(K)`:

    ℓ^{e,i}(M)   = ∫ₑ (t⋅M⋅n) μᵢ ds,     i = 0 … r
    ℓ^{K,e,q}(M) = ∫_K (t_e⋅M⋅n_e) q dK

"""
function GLSRefFE(::Type{T}, p::Polytope{D}, order::Integer) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The Gopalakrishnan-Lederer-Schoberl element is only defined on 2D simplices, got $p.
  """
  @notimplementedif order < 0 """\n
  The Gopalakrishnan-Lederer-Schoberl element needs order ≥ 0, got $order.
  """
  M = TensorValue{2,2,T,4}
  prebasis = BernsteinBasisOnSimplex(Val(2), M, order)

  # DoF moments
  fb = LegendreBasis(Val(1), T, order)                       # μ₀..μ_r on an edge
  cb = order > 0 ? BernsteinBasisOnSimplex(Val(2), T, order - 1) : nothing
  ntmom(φ, μ, ds) = Broadcasting(Operation(*))(              # ∫ₑ (t⋅M⋅n) μ ds
    Broadcasting(Operation(⊙))(φ, _gls_tn(ds)), μ
  )
  # the interior moments test with each edge's own (t⊗n), a fixed reference
  # tensor, over the cell; `E` below is that tensor
  cmom(E) = (φ, μ, ds) -> Broadcasting(Operation(*))(
    Broadcasting(Operation(⊙))(φ, ConstantField(E)), μ
  )

  # Constraint moments: tr(M) against a complete P_r(K)
  qb = BernsteinBasisOnSimplex(Val(2), T, order)
  trmom(φ, μ, ds) = Broadcasting(Operation(*))(Broadcasting(Operation(tr))(φ), μ)

  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  ts, ns = _edge_frames(p)
  moments = Tuple[(edges, ntmom, fb)]
  if order > 0
    append!(moments, [(cell, cmom(outer(ts[e], ns[e])), cb) for e in 1:num_faces(p, 1)])
  end
  push!(moments, (cell, trmom, qb))                          # the constraints

  full = MomentBasedDofBasis(p, prebasis, moments)
  ndofs = 3 * (order + 1) * (order + 2) ÷ 2

  @check length(full) == length(prebasis) """\n
  The augmented element needs as many functionals as prebasis functions, got
  $(length(full)) and $(length(prebasis)).
  """
  # the DoFs come first -- edges are enumerated before the cell, and within the
  # cell the interior moments before the constraints -- so the slice suffices
  shapefuns = linear_combination(inv(evaluate(full, prebasis))[:, 1:ndofs], prebasis)
  dofs = restrict(full, collect(1:ndofs))

  @check maximum(abs, evaluate(dofs, shapefuns) - Matrix{Float64}(I, ndofs, ndofs)) < 1e-10 """\n
  The augmented GLS Vandermonde is singular or badly conditioned.
  """
  GenericRefFE{GopalakrishnanLedererSchoberl}(
    ndofs, p, prebasis, dofs, DivConformity(), nothing,
    get_face_own_moments(dofs), shapefuns
  )
end

# The tensor t ⊗ n of the current edge of `ds`, so that M ⊙ (t⊗n) = t⋅M⋅n.
_gls_tn(ds) = ConstantField(outer(get_edge_tangent(ds).value,
                                  _rot90(get_edge_tangent(ds).value)))

function ReferenceFE(p::Polytope, ::GopalakrishnanLedererSchoberl, ::Type{T}, order) where T
  GLSRefFE(T, p, order)
end

# edge DoFs only flip sign under reversal, cell DoFs are permutation-invariant
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{GopalakrishnanLedererSchoberl}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end
