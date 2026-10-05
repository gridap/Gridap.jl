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
    ArnoldWintherCRefFE(::Type{T}, K::Polytope{2})

The conforming Arnold--Winther reference FE on the triangle `K`, with `T` the
scalar type: 24 DoFs and, writing `S` for the symmetric 2×2 tensors,

    AWc(K) = {τ ∈ P₃(K;S) : div τ ∈ P₁(K;R²)}

of dimension 24 [Arnold & Winther, Numer. Math. 92 (2002) 401].

This element is divergence conforming in the sense that the normal-normal trace
on facets `n⋅τ⋅n` is pointwise continuous.

The implementation follows the augmented element approach of [Kirby, SMAI-JCM 4 (2018) 197].

# Extended help

## Prebasis

The prebasis is taken as `P₃(K;S)`, of dimension 30, with constraints

    ∫_K (div τ)_c q dK = 0    ∀ q ∈ P₂(K) ∩ P₁(K)^⊥,  c = 1, 2

enforced using moments, yielding the 24 DoFs of the conforming Arnold--Winther element.

## Moments

At each vertex `v` the three independent components,

    ℓ^{v,c}(τ) = τ(v) ⊙ E_c,    E_c ∈ {e₁⊗e₁, e₁⊗e₂, e₂⊗e₂}

per edge `e`, with unit tangent `t`, normal `n = R t` and `μᵢ` the
L²(e)-orthonormal Legendre basis,

    ℓ^{nn,i}_e(τ) = ∫ₑ (n⋅τ⋅n) μᵢ ds,    ℓ^{nt,i}_e(τ) = ∫ₑ (n⋅τ⋅t) μᵢ ds,   i = 0, 1

and over the cell the same three components,

    ℓ^{K,c}(τ) = ∫_K τ ⊙ E_c dK.

"""
function ArnoldWintherCRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The conforming Arnold-Winther element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), SymTensorValue{2,T}, 3)

  # DoF moments
  fb = LegendreBasis(Val(1), T, 1) # μ₀, μ₁ on the edge
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
  Emom(φ, μ, ds) = Broadcasting(Operation(⊙))(φ, μ)  # ∫_K τ⊙μ dK

  verts = get_dimrange(p, 0)
  edges = get_dimrange(p, 1)
  cell = get_dimrange(p, 2)
  moments = Tuple[
    (verts, Emom, Ei),                      # Vertex moments
    (edges, nnmom, fb), (edges, ntmom, fb), # Edge moments
    (cell, Emom, Ei),                       # Cell moments
  ]

  # Constraint moments. qb is a basis for P₂(K) ∩ P₁(K)^⊥ using hierarchical
  # Dubiner, for all components
  qb = DubinerBasis(Val(2), VectorValue{2,T}, 2, Polynomials._p_complement_filter(1))
  divmom(divφ, μ, ds) = Broadcasting(Operation(⊙))(divφ, μ)  # ∫_K (div τ)⊙μ dK
  constraints = Tuple[(cell, divmom, qb)]

  dofs = MomentBasedDofBasis(p, prebasis, moments)
  cons = MomentBasedDofBasis(p, prebasis, constraints, divergence)
  full = vcat(dofs, cons)
  ndofs = length(dofs)

  @check length(full) == length(prebasis) """\n
  The augmented element needs as many functionals as prebasis functions, got
  $(length(full)) and $(length(prebasis)).
  """
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

# vertex/cell DoFs are permutation-invariant, edge DoFs only flip sign under reversal
function get_face_own_dofs_permutations(
  reffe::GenericRefFE{ArnoldWintherC}, conf::Conformity
)
  _identity_dof_permutations(reffe, conf)
end

################################################################################
# Change of basis
#
# The cell-local map. The mesh-level `compute_cell_bases_changes`, which reads
# the orientation data off the model and maps this over the cells, lives in
# src/FESpaces/Pullbacks.jl.

#     AWCChangeOfBasis(reffe, transposed_inverse)
#
# A 3×3 block per vertex, a 4×4 block per edge, and the identity on the interior.
struct AWCChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  vertex_dofs::Vector{Vector{Int}}
  edge_dofs::Vector{Vector{Int}}
  ndofs::Int
  transposed_inverse::Bool
end

function AWCChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  ts, ns = _edge_frames(p)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  vertex_dofs = [own[v] for v in 1:nv]
  edge_dofs = [own[nv+e] for e in 1:num_faces(p, 1)]
  AWCChangeOfBasis(ts, ns, vertex_dofs, edge_dofs, num_dofs(reffe), transposed_inverse)
end

function return_cache(k::AWCChangeOfBasis, Jt, σ)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::AWCChangeOfBasis, Jt, σ)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0     # the interior DoFs are left as the push-forward
  end

  J = transpose(Jt)
  detJ = det(Jt)

  # W carries det(J)⁻² _congruence_matrix(J) on the vertices, so W⁻¹ carries
  # det(J)² _congruence_matrix(J⁻¹) and Wᵀ the transpose of the former.
  Bv = k.transposed_inverse ? transpose(_congruence_matrix(J)) / detJ^2 :
       _congruence_matrix(inv(J)) * detJ^2
  for dofs in k.vertex_dofs
    for i in 1:3, j in 1:3
      M[dofs[i], dofs[j]] = Bv[i, j]
    end
  end

  _aw_edge_blocks!(M, k.tangents, k.normals, k.edge_dofs, Jt, σ, k.transposed_inverse)
  return M
end

################################################################################
# DOF scaling
#
function get_dofscale_setter_function(
  reffe::GenericRefFE{ArnoldWintherC}, ::DoubleContraVariantPiolaMap
)
  _face_dim_dofscale_setter(reffe)
end

