"""
    struct Argyris <: ReferenceFEName end

Reference FE name for the quintic Argyris triangle, a C¹-conforming element. See
[`ArgyrisRefFE`](@ref). Its singleton instance is [`argyris`](@ref).
"""
struct Argyris <: ReferenceFEName end

"""
    const argyris = Argyris()

Singleton of the [`Argyris`](@ref) reference FE name.
"""
const argyris = Argyris()

Pushforward(::Type{Argyris}) = IdentityPiolaMap()

"""
    ArgyrisRefFE(::Type{T}, K::Polytope{2})

The Argyris reference FE on the triangle `K`, with `T` the scalar type: the
quintic C¹ plate element with 21 DoFs of [Argyris, Fried & Scharpf, Aero. J. 72 (1968) 701].

The implementation follows [Kirby, SMAI-JCM 4 (2018) 197, §3.3.2].

# Extended help

## Prebasis

The prebasis is taken as `P₅(K)`, of dimension 21.

## Moments

At each vertex `v` the value, the two gradient components and the three
independent Hessian components; per edge `e` with unit tangent `t` and normal
`n = R t`, the mean normal derivative:

    ℓ^{v,0}(u)   = u(v)
    ℓ^{v,i}(u)   = ∇u(v) ⋅ eᵢ,               i = 1, 2
    ℓ^{v,ij}(u)  = ∇∇u(v) ⊙ E_ij,            E_ij ∈ {e₁⊗e₁, e₁⊗e₂, e₂⊗e₂}
    ℓ^e(u)       = ∫ₑ (∇u⋅n) ds

"""
function ArgyrisRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The Argyris element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), T, 5)  # P₅(K), dim 21
  dofs = _argyris_dof_basis(T, p, prebasis)
  face_own_dofs = _argyris_face_own_dofs(p)

  ndofs = 6*num_faces(p, 0) + num_faces(p, 1)
  @check length(dofs) == ndofs == length(prebasis)

  GenericRefFE{Argyris}(ndofs, p, prebasis, dofs, H2Conformity(), nothing, face_own_dofs)
end

function ReferenceFE(p::Polytope, ::Argyris, ::Type{T}) where T
  ArgyrisRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::Argyris, ::Type{T}, order) where T
  @notimplementedif order != 5 """\n
  The Argyris element exists for order 5 only, got $order.
  """
  ArgyrisRefFE(T, p)
end

function _argyris_dof_basis(::Type{T}, p::Polytope{2}, prebasis) where T
  verts = get_dimrange(p, 0)
  edges = get_dimrange(p, 1)

  vmom(φ, μ, ds) = Broadcasting(Operation(*))(φ, μ)
  c0 = MonomialBasis(Val(0), T, 0)  # the constant on a vertex
  values = MomentBasedDofBasis(p, prebasis, Tuple[(verts, vmom, c0)])

  c0_grad = map(
    constant_field, representatives_of_componentbasis_dual(VectorValue{2,T})
  )
  c0_hess = map(
    constant_field, representatives_of_componentbasis_dual(SymTensorValue{2,T})
  )
  Dmom(Dφ, μ, ds) = Broadcasting(Operation(⊙))(Dφ, μ)  # ∇u(v)⊙μ or ∇∇u(v)⊙μ
  grads = MomentBasedDofBasis(p, prebasis, Tuple[(verts, Dmom, c0_grad)], ∇)
  hessians = MomentBasedDofBasis(p, prebasis, Tuple[(verts, Dmom, c0_hess)], ∇∇)

  e0 = MonomialBasis(Val(1), T, 0)  # the constant on an edge
  function emom(φ, μ, ds)  # σ_e(∇u,μ) = ∫ₑ (∇u⋅n) μ ds
    n = _edge_normal(ds)
    Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φ, n), μ)
  end
  normals = MomentBasedDofBasis(p, prebasis, Tuple[(edges, emom, e0)], ∇)

  return vcat(values, grads, hessians, normals)
end

function _argyris_face_own_dofs(p::Polytope{2})
  nv, ne = num_faces(p, 0), num_faces(p, 1)
  own = [Int[] for _ in 1:num_faces(p)]
  for v in 1:nv
    own[v] = [v, nv + 2*v - 1, nv + 2*v, 3*nv + 3*v - 2, 3*nv + 3*v - 1, 3*nv + 3*v]
  end
  for e in 1:ne
    own[nv+e] = [6*nv + e]
  end
  return own
end

function get_face_own_dofs_permutations(reffe::GenericRefFE{Argyris}, conf::Conformity)
  _identity_dof_permutations(reffe, conf)
end

################################################################################
# Change of basis

#     _congruence_matrix(A) -> TensorValue{3,3}
#
# The 3×3 tensor of the congruence `H ↦ A H Aᵀ` acting on symmetric 2×2 tensor,
# in the coordinates `(H₁₁, H₁₂, H₂₂)` — equivalently the second symmetric power
# `Sym²(A)`, i.e. `A ⊗ A` restricted to the symmetric subspace.
function _congruence_matrix(A)
  a11, a12, a21, a22 = A[1,1], A[1,2], A[2,1], A[2,2]
  # column-major, one column per basis tensor: [1 0;0 0], [0 1;1 0], [0 0;0 1]
  return TensorValue{3,3}(
    a11*a11, a11*a21,             a21*a21,
    2*a11*a12, a11*a22 + a12*a21, 2*a21*a22,
    a12*a12, a12*a22,             a22*a22
  )
end

# Argyris is mapped by the plain pullback u = û∘F⁻¹. Writing K = J⁻ᵀ, and using
# that F is affine on a simplex so no second derivative of F appears,
#
#   ∇u = K ∇û,        D²u = K D²û Kᵀ,
#
# the vertex value DoFs are preserved, the vertex gradient DoFs mix within a
# vertex through K, and the vertex Hessian DoFs mix within a vertex through
# `_congruence_matrix(K)`. The edge DoFs behave as in Morley: with G = (JᵀJ)⁻¹
# and n = R t,
#
#   F∗(δₑᵏ) = Aₖ δ̂ₑᵏ + Bₖ (δ̂ᵥᵇ - δ̂ᵥᵃ),   Aₖ = det(J) n̂ᵀGn̂,  Bₖ = det(J) t̂ᵀGn̂,
#
# coupling each edge row to the two *value* DoFs of its endpoints. So
#
#         ⎡ I    0    0    0 ⎤                  ⎡ I       0     0      0   ⎤
#   W  =  ⎢ 0    Kg   0    0 ⎥ ,      W⁻¹  =    ⎢ 0       Kg⁻¹  0      0   ⎥
#         ⎢ 0    0    Kh   0 ⎥                  ⎢ 0       0     Kh⁻¹   0   ⎥
#         ⎣ Bv   0    0    A ⎦                  ⎣ -A⁻¹Bv  0     0      A⁻¹ ⎦
#
# with Kg, Kh block diagonal over the vertices and A = diag(Aₖ) — block
# triangular, not block diagonal. Kg⁻¹ and Kh⁻¹ are the same constructions
# applied to K⁻¹ = Jᵀ, so no tensor is ever inverted numerically.
#
# Edge orientation follows `_edge_signs`, with D = diag(1,…,1,σ₁,σ₂,σ₃) folded in
# as P = W⁻¹D and P⁻ᵀ = WᵀD. The vertex DoFs need no such treatment, being stated
# in the global Cartesian frame and so the same functional for every cell that
# touches the vertex.

#     ArgyrisChangeOfBasis(p, transposed_inverse)
#
# Builds, from the (transposed) Jacobian of a cell's geometrical map and that
# cell's edge orientation signs `σ`, either the change of basis `P`
# (`transposed_inverse = false`) or `P⁻ᵀ` (`transposed_inverse = true`).
struct ArgyrisChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  edge_vertices::Vector{Vector{Int}}
  transposed_inverse::Bool
end

function ArgyrisChangeOfBasis(p::Polytope{2}, transposed_inverse::Bool)
  ts, ns = _edge_frames(p)
  ArgyrisChangeOfBasis(ts, ns, get_faces(p, 1, 0), transposed_inverse)
end

function return_cache(k::ArgyrisChangeOfBasis, Jt, σ)
  nv, ne = length(k.tangents), length(k.edge_vertices)
  CachedArray(zeros(Float64, 6*nv + ne, 6*nv + ne))
end

function evaluate!(cache, k::ArgyrisChangeOfBasis, Jt, σ)
  nv = length(k.tangents)   # a triangle: as many vertices as edges
  ne = length(k.edge_vertices)
  ndofs = 6*nv + ne
  setsize!(cache, (ndofs, ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))

  detJ = det(Jt)
  G = inv(Jt ⋅ transpose(Jt))   # (JᵀJ)⁻¹, since Jt = Jᵀ

  # W carries K = J⁻ᵀ on the gradients and `_congruence_matrix(K)` on the
  # Hessians, so W⁻¹ carries K⁻¹ = Jᵀ = Jt and `_congruence_matrix(Jt)`, while Wᵀ
  # carries their transposes.
  K = inv(Jt)
  Kg = ifelse(k.transposed_inverse, transpose(K), Jt)
  Kh = k.transposed_inverse ? transpose(_congruence_matrix(K)) :
       _congruence_matrix(Jt)

  gof = nv           # offset of the gradient DoFs
  hof = 3*nv         # offset of the Hessian DoFs
  eof = 6*nv         # offset of the edge DoFs

  for v in 1:nv
    M[v, v] = 1.0
    for i in 1:2, j in 1:2
      M[gof+2*v-2+i, gof+2*v-2+j] = Kg[i, j]
    end
    for i in 1:3, j in 1:3
      M[hof+3*v-3+i, hof+3*v-3+j] = Kh[i, j]
    end
  end

  for e in 1:ne
    t̂, n̂ = k.tangents[e], k.normals[e]
    Gn̂ = G ⋅ n̂
    A = detJ * (n̂ ⋅ Gn̂)
    B = detJ * (t̂ ⋅ Gn̂)
    σe = σ[e]
    va, vb = k.edge_vertices[e]

    # D = diag(1,…,1,σ…) scales the last columns of the block below.
    if k.transposed_inverse
      # WᵀD
      M[eof+e, eof+e] = A * σe
      M[va, eof+e] = -B * σe
      M[vb, eof+e] = B * σe
    else
      # W⁻¹D
      M[eof+e, eof+e] = σe / A
      M[eof+e, va] = B / A
      M[eof+e, vb] = -B / A
    end
  end

  return M
end
