"""
    struct HellanHerrmannJohnson <: ReferenceFEName end

Reference FE name for the Hellan--Herrmann--Johnson element, a
normal-normal-continuous symmetric-tensor-valued element for the Kirchhoff
plate. See [`HHJRefFE`](@ref). Its singleton instance is [`hhj`](@ref).
"""
struct HellanHerrmannJohnson <: ReferenceFEName end

"""
    const hhj = HellanHerrmannJohnson()

Singleton of the [`HellanHerrmannJohnson`](@ref) reference FE name.
"""
const hhj = HellanHerrmannJohnson()

Pushforward(::Type{HellanHerrmannJohnson}) = DoubleContraVariantPiolaMap()

"""
    HHJRefFE(::Type{T}, K::Polytope{2}, order::Integer)

The Hellan--Herrmann--Johnson reference FE of degree `r = order` on the triangle
`K`, with `T` the scalar type [Arnold & Walker, SIAM J. Numer. Anal. 58 (2020)
2829, (2.3)]. Writing `S` for the symmetric 2×2 matrices,

    HHJ_r(K) = P_r(K;S)

of dimension `3(r+1)(r+2)/2`, for any `r ≥ 0`.

# Extended help

## Prebasis

The prebasis is taken as `P_r(K;S)`, of dimension `3(r+1)(r+2)/2`.

## Moments

Per edge `e` with unit normal `n` and `μᵢ` the L²(e)-orthonormal Legendre basis
of `P_r(e)`, the normal--normal moments; over the cell, the moments against a
basis `{τ}` of `P_{r-1}(K;S)`:

    ℓ^{e,i}(ϕ) = ∫ₑ (n⋅ϕ⋅n) μᵢ ds,    i = 0 … r
    ℓ^{K,τ}(ϕ) = ∫_K ϕ ⊙ τ dK.

"""
function HHJRefFE(::Type{T}, p::Polytope{D}, order::Integer) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The Hellan-Herrmann-Johnson element is only defined on 2D simplices, got $p.
  """
  @notimplementedif order < 0 """\n
  The Hellan-Herrmann-Johnson element needs order ≥ 0, got $order.
  """
  S = SymTensorValue{2,T}

  prebasis = BernsteinBasisOnSimplex(Val(2), S, order)
  fb = LegendreBasis(Val(1), T, order)  # μ₀..μ_r on the edge, μᵢ of parity (-1)ⁱ
  cb = order > 0 ? BernsteinBasisOnSimplex(Val(2), S, order - 1) : nothing

  function fmom(φ, μ, ds)  # σ_e(ϕ,μ) = ∫ₑ (n⋅ϕ⋅n) μ ds
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ, n)
    φnn = Broadcasting(Operation(⋅))(φn, n)
    Broadcasting(Operation(*))(φnn, μ)
  end
  cmom(φ, μ, ds) = Broadcasting(Operation(⊙))(φ, μ)  # σ_K(ϕ,τ) = ∫_K ϕ⊙τ dK

  moments = Tuple[(get_dimrange(p, 1), fmom, fb)]
  if order > 0
    push!(moments, (get_dimrange(p, 2), cmom, cb))
  end

  MomentBasedReferenceFE(hhj, p, prebasis, moments, DivConformity())
end

function ReferenceFE(p::Polytope, ::HellanHerrmannJohnson, ::Type{T}, order) where T
  HHJRefFE(T, p, order)
end

# edge DoFs only flip sign under reversal, cell DoFs are permutation-invariant
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{HellanHerrmannJohnson}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end
