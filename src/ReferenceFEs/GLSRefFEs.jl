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

# The tensor t ⊗ n of the current edge of `ds`, so that M ⊙ (t⊗n) = t⋅Mn. Here
# n = R t, which is what makes the pair flip together under a reversal of the
# edge and the functional invariant.
_gls_tn(ds) = ConstantField(outer(get_edge_tangent(ds).value,
                                  _rot90(get_edge_tangent(ds).value)))

"""
    GLSRefFE(::Type{T}, p::Polytope{2}, order::Integer)

The Gopalakrishnan--Lederer--Schoberl reference FE of the second kind, degree
`r = order` on the triangle `K`, with `T` the scalar type. Writing `M₀` for the
traceless 2×2 matrices,

    GLS_r(K) = P_r(K;M₀)

of dimension `3(r+1)(r+2)/2`, for any `r ≥ 0`
[Gopalakrishnan, Lederer & Schoberl, SIAM J. Numer. Anal. 58 (2020) 706].

Assembled as the augmented element of [Kirby, SMAI-JCM 4 (2018) 197, §5.5]: the
DoFs *and* the constraint functionals are declared together over an ambient
prebasis they form a dual basis for, the generalized Vandermonde is inverted, and
the columns dual to the DoFs are kept — they lie in the constrained space,
being annihilated by every constraint. Note `length(get_prebasis(reffe))` is
therefore larger than `num_dofs(reffe)`.

## Prebasis

`P_r(K;M)`, the *full* matrix-valued space of dimension `4(r+1)(r+2)/2` — the
ambient space, not `GLS_r(K)` itself.

## Moments

Per edge `e` with unit tangent `t`, normal `n = R t` and `μᵢ` the
L²(e)-orthonormal Legendre basis of `P_r(e)`; and over the cell, the same
normal-tangential functional of each edge's frame against a basis `{q}` of
`P_{r-1}(K)`:

    ℓ^{e,i}(M)   = ∫ₑ (t⋅Mn) μᵢ ds,     i = 0 … r
    ℓ^{K,e,q}(M) = ∫_K (t_e⋅Mn_e) q dK

`3(r+1) + 3r(r+1)/2 = 3(r+1)(r+2)/2`. The edge DoFs are the ones that glue: they
give the normal--tangential continuity the MCS method needs, while the full
matrix jump stays O(1). At `r = 0` there are no interior moments.

`t⋅Mn = M ⊙ (t⊗n)` is bilinear with both directions flipping under a reversal of
the edge, so it is invariant and no orientation convention is needed.

## Constraints

    ∫_K tr(M) q dK = 0    ∀ q ∈ P_r(K)      ((r+1)(r+2)/2 of them)

`tr M` lies in `P_r(K)` a priori, and the moment runs over a *complete* `P_r(K)`,
so this says `tr M ≡ 0`. `3(r+1)(r+2)/2 + (r+1)(r+2)/2 = 4(r+1)(r+2)/2`.
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
  ntmom(φ, μ, ds) = Broadcasting(Operation(*))(              # ∫ₑ (t⋅Mn) μ ds
    Broadcasting(Operation(⊙))(φ, _gls_tn(ds)), μ
  )
  # the interior moments test with each edge's own (t⊗n), a fixed reference
  # tensor, over the cell; `E` below is that tensor
  cmom(E) = (φ, μ, ds) -> Broadcasting(Operation(*))(
    Broadcasting(Operation(⊙))(φ, ConstantField(E)), μ
  )

  # Constraint moments: tr(M) = M ⊙ I against a complete P_r(K)
  qb = BernsteinBasisOnSimplex(Val(2), T, order)
  Id = TensorValue(one(T), zero(T), zero(T), one(T))

  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  ts, ns = _edge_frames(p)
  moments = Tuple[(edges, ntmom, fb)]
  if order > 0
    append!(moments, [(cell, cmom(outer(ts[e], ns[e])), cb) for e in 1:num_faces(p, 1)])
  end
  push!(moments, (cell, cmom(Id), qb))                       # the constraints

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

function ReferenceFE(p::Polytope, ::GopalakrishnanLedererSchoberl, ::Type{T}, order) where T
  GLSRefFE(T, p, order)
end

# Identity for every admissible vertex permutation of every face: the DoFs of an
# edge are ordered by the degree of their Legendre weight, which both adjacent
# cells agree on, and reversing the edge only changes the signs of the odd ones,
# which the change of basis carries.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{GopalakrishnanLedererSchoberl}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end
