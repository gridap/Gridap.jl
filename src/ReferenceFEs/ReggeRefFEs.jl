"""
    struct Regge <: ReferenceFEName end

Reference FE name for the Regge element, a tangential-tangential-continuous
symmetric-tensor-valued element. See [`ReggeRefFE`](@ref). Its singleton instance
is [`regge`](@ref).
"""
struct Regge <: ReferenceFEName end

"""
    const regge = Regge()

Singleton of the [`Regge`](@ref) reference FE name.
"""
const regge = Regge()

Pushforward(::Type{Regge}) = DoubleCoVariantPiolaMap()

"""
    ReggeRefFE(::Type{T}, K::Polytope{2}, order::Integer)

The Regge reference FE of degree `r = order` on the triangle `K`, with `T` the
scalar type. Writing `S` for the symmetric 2×2 tensors,

    Regge_r(K) = P_r(K;S)

of dimension `3(r+1)(r+2)/2`, for any `r ≥ 0`. This element is H(curl-curl)
conform, since the tangential-tangential trace `t⋅M⋅t` is continuous on the
boundary.

# Extended help

## Prebasis

`P_r(K;S)`, of dimension `3(r+1)(r+2)/2`, the element's space itself.

## Moments

Per edge `e` with unit tangent `t` and `μᵢ` the L²(e)-orthonormal Legendre basis
of `P_r(e)`, the tangential--tangential moments; over the cell, the moments
against a basis `{τ}` of `P_{r-1}(K;S)`:

    ℓ^{e,i}(M) = ∫ₑ (t⋅M⋅t) μᵢ ds,    i = 0 … r
    ℓ^{K,τ}(M) = ∫_K M ⊙ τ dK

`3(r+1) + 3r(r+1)/2 = 3(r+1)(r+2)/2`. The edge DoFs are the ones that glue: they
give the tangential--tangential continuity Regge needs, while the full tensor
jump stays O(1). At `r = 0` there are no interior moments.

`t⋅M⋅t` is quadratic in `t`, so its sign is irrelevant and no orientation
convention is needed at all — the same freedom [`HHJRefFE`](@ref) has with
`n⋅ϕ⋅n`.
"""
function ReggeRefFE(::Type{T}, p::Polytope{D}, order::Integer) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The Regge element is only defined on 2D simplices, got $p.
  """
  @notimplementedif order < 0 """\n
  The Regge element needs order ≥ 0, got $order.
  """
  S = SymTensorValue{2,T}

  prebasis = BernsteinBasisOnSimplex(Val(2), S, order)
  fb = LegendreBasis(Val(1), T, order)  # μ₀..μ_r on the edge, μᵢ of parity (-1)ⁱ
  cb = order > 0 ? BernsteinBasisOnSimplex(Val(2), S, order - 1) : nothing

  function fmom(φ, μ, ds)  # σ_e(M,μ) = ∫ₑ (t⋅M⋅t) μ ds
    t = get_edge_tangent(ds)
    Mt = Broadcasting(Operation(⋅))(φ, t)
    Mtt = Broadcasting(Operation(⋅))(Mt, t)
    Broadcasting(Operation(*))(Mtt, μ)
  end
  cmom(φ, μ, ds) = Broadcasting(Operation(⊙))(φ, μ)  # σ_K(M,τ) = ∫_K M⊙τ dK

  moments = Tuple[(get_dimrange(p, 1), fmom, fb)]
  if order > 0
    push!(moments, (get_dimrange(p, 2), cmom, cb))
  end

  MomentBasedReferenceFE(regge, p, prebasis, moments, DivConformity())
end

function ReferenceFE(p::Polytope, ::Regge, ::Type{T}, order) where T
  ReggeRefFE(T, p, order)
end

# edge DoFs only flip sign under reversal, cell DoFs are permutation-invariant
function get_face_own_dofs_permutations(reffe::GenericRefFE{Regge}, conf::Conformity)
  _identity_dof_permutations(reffe, conf)
end

################################################################################
# DOF scaling
#
# The change of basis is `EdgeScalingChangeOfBasis`, in Pullbacks.jl.
function get_dofscale_setter_function(reffe::GenericRefFE{Regge}, ::DoubleCoVariantPiolaMap)
  _face_dim_dofscale_setter(reffe)
end

