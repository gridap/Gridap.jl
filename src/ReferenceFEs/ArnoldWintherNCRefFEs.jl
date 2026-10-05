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
Pushforward(::Type{ArnoldWintherNC}, ::L2Conformity) = DoubleContraVariantPiolaMap()

"""
    ArnoldWintherNCRefFE(::Type{T}, K::Polytope{2})

The nonconforming Arnold--Winther reference FE on the triangle `K`, with `T` the
scalar type: 15 DoFs and, writing `S` for the symmetric 2×2 tensors,

    AWnc(K) = {τ ∈ P₂(K;S) : (n⋅τ⋅n)|ₑ ∈ P₁(e) ∀ e ⊂ ∂K}

of dimension 15 [Arnold & Winther, M3AS 13 (2003) 295].

This element is not divergence conforming so it's conformity is `:L2`. But the
normal-normal moment on facets, `∫_f n⋅τ⋅n df`, are preserved, that is single
valued on both sides.

The implementation follows the augmented element approach of [Kirby, SMAI-JCM 4 (2018) 197].

# Extended help

## Prebasis

The prebasis is taken as `P₂(K;S)`, of dimension 18, with constraints

    ∫ₑ (n⋅τ⋅n) μ ds = 0    ∀ μ ∈ P₂(e) ∩ P₁(e)^⊥,  ∀ e ⊂ ∂K

enforced using moments, yielding the 15 DoFs of the nonconforming Arnold--Winther element.

## Moments

Per edge `e`, with unit tangent `t`, normal `n = R t` and `μᵢ` the
L²(e)-orthonormal Legendre basis:

    ℓ^{nn,i}_e(τ) = ∫ₑ (n⋅τ⋅n) μᵢ ds,    i = 0, 1
    ℓ^{nt,i}_e(τ) = ∫ₑ (n⋅τ⋅t) μᵢ ds,    i = 0, 1

four per edge; and over the cell the three independent components,

    ℓ^{K,c}(τ) = ∫_K τ ⊙ E_c dK,    E_c ∈ {e₁⊗e₁, e₁⊗e₂, e₂⊗e₂}.

"""
function ArnoldWintherNCRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The nonconforming Arnold-Winther element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), SymTensorValue{2,T}, 2)
  fb = LegendreBasis(Val(1), T, 1)                            # μ₀, μ₁ : the DoF weights
  # μ₂ : the constraint, the degree-2 Legendre polynomial alone. It spans
  # P₂(e) ∩ P₁(e)^⊥, so the single moment against it states (n⋅τ⋅n)|ₑ ∈ P₁(e).
  gb = LegendreBasis(Val(1), T, 2, Polynomials._p_complement_filter(1))
  Ei = map(constant_field, representatives_of_componentbasis_dual(SymTensorValue{2,T}))

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
  cmom(φ, μ, ds) = Broadcasting(Operation(⊙))(φ, μ)  # ∫_K τ⊙μ dK

  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  moments = Tuple[
    (edges, nnmom, fb), (edges, ntmom, fb), (edges, nnmom, gb), # Edge moments
    (cell, cmom, Ei),                                           # Cell moments
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
    ndofs, p, prebasis, dofs, L2Conformity(), nothing,
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

# edge DoFs only flip sign under reversal, cell DoFs are permutation-invariant
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{ArnoldWintherNC}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end

################################################################################
# Change of basis
#
# The cell-local map. The mesh-level `compute_cell_bases_changes`, which reads
# the orientation data off the model and maps this over the cells, lives in
# src/FESpaces/Pullbacks.jl.

# Both AW elements are mapped by the double contravariant Piola map
# τ = det(J)⁻² J τ̂ Jᵀ, and both carry the same four edge DoFs -- the degree 0 and
# 1 moments of n⋅τ⋅n and of n⋅τ⋅t. Those give one 4×4 block per edge, in the DoF
# order (nn0, nn1, nt0, nt1),
#
#   W = [1/L  0    0  0;  0  1/L  0  0;  α  0  β  0;  0  α  0  β],
#
# with L = ‖J t̂ₑ‖, α = ã/(det(J) L), β = L/det(J) and ã = n̂ᵀ(JᵀJ)t̂. The interior
# DoFs are left as the push-forward of the reference ones -- they are cell-owned
# and shared with nobody -- which sidesteps a dense interior block entirely, as
# Gridap already does for the cell moments of Raviart-Thomas.
#
# The conforming element adds a 3×3 block per vertex, det(J)⁻² times the tensor
# of H ↦ J H Jᵀ, i.e. `_congruence_matrix` with A = J rather than the A = J⁻ᵀ the
# Argyris Hessian block uses. Both DoF kinds are invariant under reversing an
# edge -- n⋅τ⋅n is quadratic in n, n⋅τ⋅t bilinear with both flipping -- so only the
# parity of the Legendre weight enters σ.

function _aw_edge_blocks!(M, tangents, normals, edge_dofs, Jt, σ, transposed_inverse)
  detJ = det(Jt)
  JtJ = Jt ⋅ transpose(Jt)   # JᵀJ, since Jt = Jᵀ

  for e in eachindex(edge_dofs)
    t̂, n̂ = tangents[e], normals[e]
    L = norm(t̂ ⋅ Jt)               # ‖J t̂‖
    α = (n̂ ⋅ (JtJ ⋅ t̂)) / (detJ * L)
    β = L / detJ
    reversed = σ[e] < 0

    # DoF order within an edge: (nn,0), (nn,1), (nt,0), (nt,1)
    dofs = edge_dofs[e]
    nmom = length(dofs) ÷ 2
    for i in 1:nmom
      s = ifelse(reversed && isodd(i - 1), -1.0, 1.0)   # parity of the Legendre weight
      dnn, dnt = dofs[i], dofs[nmom+i]
      if transposed_inverse
        # WᵀD
        M[dnn, dnn] = s / L
        M[dnn, dnt] = s * α
        M[dnt, dnt] = s * β
      else
        # W⁻¹D
        M[dnn, dnn] = s * L
        M[dnt, dnn] = -s * α * L / β
        M[dnt, dnt] = s / β
      end
    end
  end
end

#     AWNCChangeOfBasis(reffe, transposed_inverse)
#
# One 4×4 block per edge, and the identity on the interior.
struct AWNCChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  edge_dofs::Vector{Vector{Int}}
  ndofs::Int
  transposed_inverse::Bool
end

function AWNCChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  ts, ns = _edge_frames(p)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  edge_dofs = [own[nv+e] for e in 1:num_faces(p, 1)]
  AWNCChangeOfBasis(ts, ns, edge_dofs, num_dofs(reffe), transposed_inverse)
end

function return_cache(k::AWNCChangeOfBasis, Jt, σ)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::AWNCChangeOfBasis, Jt, σ)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0     # the interior DoFs are left as the push-forward
  end
  _aw_edge_blocks!(M, k.tangents, k.normals, k.edge_dofs, Jt, σ, k.transposed_inverse)
  return M
end

################################################################################
# DOF scaling
#
function get_dofscale_setter_function(
  reffe::GenericRefFE{ArnoldWintherNC}, ::DoubleContraVariantPiolaMap
)
  _face_dim_dofscale_setter(reffe)
end

