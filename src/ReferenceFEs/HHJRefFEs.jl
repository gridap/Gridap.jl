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
    HHJRefFE(::Type{T}, p::Polytope{2}, order::Integer)

The Hellan--Herrmann--Johnson reference FE of degree `r = order` on the triangle
`K`, with `T` the scalar type [Arnold & Walker, SIAM J. Numer. Anal. 58 (2020)
2829, (2.3)]. Writing `S` for the symmetric 2×2 matrices,

    HHJ_r(K) = P_r(K;S)

of dimension `3(r+1)(r+2)/2`, for any `r ≥ 0`. Unconstrained, like `MorleyRefFE`
and `ArgyrisRefFE`: what is non-standard is that the DoFs are not preserved by
the double contravariant Piola map. They are preserved up to one positive scalar
per edge, which makes the change of basis diagonal.

## Prebasis

`P_r(K;S)`, of dimension `3(r+1)(r+2)/2`, the element's space itself.

## Moments

Per edge `e` with unit normal `n` and `μᵢ` the L²(e)-orthonormal Legendre basis
of `P_r(e)`, the normal--normal moments; over the cell, the moments against a
basis `{τ}` of `P_{r-1}(K;S)`:

    ℓ^{e,i}(ϕ) = ∫ₑ (n⋅ϕn) μᵢ ds,    i = 0 … r
    ℓ^{K,τ}(ϕ) = ∫_K ϕ ⊙ τ dK

`3(r+1) + 3r(r+1)/2 = 3(r+1)(r+2)/2`. The edge DoFs are the ones that glue: they
give the normal--normal continuity the HHJ method needs, while the full tensor
jump stays O(1). At `r = 0` there are no interior moments.

The normal here is the polytope's outward normal, not the `n = R t` of the
elements that need an edge frame: `n⋅ϕn` is quadratic in `n`, so its sign is
irrelevant and no orientation convention is needed at all. The one orientation
effect left is the weight — reversing an edge sends `μᵢ` to `(-1)ⁱ μᵢ` — and the
change of basis carries it.
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

  function fmom(φ, μ, ds)  # σ_e(ϕ,μ) = ∫ₑ (n⋅ϕn) μ ds
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

# Identity for every admissible vertex permutation of every face: the DoFs of an
# edge are ordered by the degree of their Legendre weight, which both adjacent
# cells agree on, and reversing the edge only changes the signs of the odd ones,
# which the change of basis carries.
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{HellanHerrmannJohnson}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end
