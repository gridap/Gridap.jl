"""
    struct Morley <: ReferenceFEName end

Reference FE name for the quadratic Morley triangle, a nonconforming element for
fourth-order problems. See [`MorleyRefFE`](@ref). Its singleton instance is
[`morley`](@ref).
"""
struct Morley <: ReferenceFEName end

"""
    const morley = Morley()

Singleton of the [`Morley`](@ref) reference FE name.
"""
const morley = Morley()

Pushforward(::Type{Morley}) = IdentityPiolaMap()

"""
    MorleyRefFE(::Type{T}, K::Polytope{2})

The Morley reference FE on the triangle `K`, with `T` the scalar type: the
quadratic nonconforming plate element of [Morley, Aero. Quart. 19 (1968) 149],
with 6 DoFs.

It's implementation conformity is `:H1`, but the element also have continuous
normal derivatives accross cell edges.

# Extended help

## Prebasis

We take prebasis `P₂(K)`, of dimension 6.

## Moments

At each vertex `v` the point value, and per edge `e` with unit tangent `t` and
normal `n = R t` (`R` the clockwise quarter turn) the mean normal derivative:

    ℓ^v(u)   = u(v)
    ℓ^e(u)   = ∫ₑ (∇u⋅n) ds

"""
function MorleyRefFE(::Type{T}, p::Polytope{D}) where {T,D}
  @notimplementedif !(D == 2 && is_simplex(p)) """\n
  The Morley element is only defined on 2D simplices, got $p.
  """
  prebasis = BernsteinBasisOnSimplex(Val(2), T, 2)  # P₂(K), dim 6
  dofs = _morley_dof_basis(T, p, prebasis)

  ndofs = num_faces(p, 0) + num_faces(p, 1)
  @check length(dofs) == ndofs == length(prebasis)

  # One DoF per vertex, then one per edge; the cell owns none
  face_own_dofs = [[i] for i in 1:ndofs]
  push!(face_own_dofs, Int[])

  GenericRefFE{Morley}(ndofs, p, prebasis, dofs, H1Conformity(), nothing, face_own_dofs)
end

function _morley_dof_basis(::Type{T}, p::Polytope{2}, prebasis) where T
  vb = MonomialBasis(Val(0), T, 0)  # the constant on a vertex
  vmom(φ, μ, ds) = Broadcasting(Operation(*))(φ, μ)      # σ_v(u,μ) = u(v) μ
  vertex_dofs = MomentBasedDofBasis(p, prebasis, Tuple[(get_dimrange(p, 0), vmom, vb)])

  eb = MonomialBasis(Val(1), T, 0)  # the constant on an edge
  function emom(φ, μ, ds)                                # σ_e(∇u,μ) = ∫(∇u⋅n)μ
    n = _edge_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ, n)
    Broadcasting(Operation(*))(φn, μ)
  end
  edge_dofs = MomentBasedDofBasis(p, prebasis, Tuple[(get_dimrange(p, 1), emom, eb)], ∇)

  return vcat(vertex_dofs, edge_dofs)
end

function ReferenceFE(p::Polytope, ::Morley, ::Type{T}) where T
  MorleyRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::Morley, ::Type{T}, order) where T
  @notimplementedif order != 2 """\n
  The Morley element exists for order 2 only, got $order.
  """
  MorleyRefFE(T, p)
end

# dofs are robust to permutation up to edge orientation, delt with in FESpaces/Pullbacks.jl
function get_face_own_dofs_permutations(reffe::GenericRefFE{Morley}, conf::Conformity)
  _identity_dof_permutations(reffe, conf)
end

################################################################################
# Change of basis
#
# The cell-local map. The mesh-level `compute_cell_bases_changes`, which reads
# the orientation data off the model and maps this over the cells, lives in
# src/FESpaces/Pullbacks.jl.

# Morley is mapped by the plain pullback u = û∘F⁻¹, under which the vertex values
# are preserved but the edge normal derivatives are not: with ∇u = J⁻ᵀ∇û,
#
#   ∇u⋅n = ∇û⋅(J⁻¹n) = a (∇û⋅n̂) + b (∇û⋅t̂),
#
# so the push-forward of an edge DoF picks up a tangential derivative, which is
# not a Morley node. With the DoFs written as moments that tangential piece is
# exactly a difference of vertex values,
#
#   ∫_{ê} ∇û⋅t̂ dŝ = û(v̂_b) - û(v̂_a),
#
# so span(N̂) is preserved after all and the transformation is a 6×6 matrix in
# closed form. With n = R t and G = (JᵀJ)⁻¹, using R J Rᵀ = det(J) J⁻ᵀ,
#
#   J⁻¹n = (det J / ‖J t̂‖) G n̂,   ds = ‖J t̂‖ dŝ,
#
# the ‖J t̂‖ cancels and, for edge k with reference endpoints v̂_a, v̂_b,
#
#   F∗(δᵥⁱ) = δ̂ᵥⁱ,
#   F∗(δₑᵏ) = Aₖ δ̂ₑᵏ + Bₖ (δ̂ᵥᵇ - δ̂ᵥᵃ),   Aₖ = det(J) n̂ᵀGn̂,  Bₖ = det(J) t̂ᵀGn̂.
#
# So W = [I 0; B A] with A = diag(Aₖ) — block *triangular*, the case Kirby notes
# for Morley and Argyris. Aₖ ≠ 0 always (n̂ᵀGn̂ > 0 by positive definiteness), and
# W⁻¹ = [I 0; -A⁻¹B A⁻¹] in closed form. Edge orientation follows `_edge_signs`,
# with D = diag(1,1,1,σ₁,σ₂,σ₃) folded in as P = W⁻¹D and P⁻ᵀ = WᵀD.

#     MorleyChangeOfBasis(p, transposed_inverse)
#
# Builds, from the (transposed) Jacobian of a cell's geometrical map and that
# cell's edge orientation signs `σ`, either the change of basis `P`
# (`transposed_inverse = false`) or `P⁻ᵀ` (`transposed_inverse = true`).
struct MorleyChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  edge_vertices::Vector{Vector{Int}}
  transposed_inverse::Bool
end

function MorleyChangeOfBasis(p::Polytope{2}, transposed_inverse::Bool)
  ts, ns = _edge_frames(p)
  MorleyChangeOfBasis(ts, ns, get_faces(p, 1, 0), transposed_inverse)
end

function return_cache(k::MorleyChangeOfBasis, Jt, σ)
  ndofs = length(k.tangents) + length(k.edge_vertices)
  CachedArray(zeros(Float64, ndofs, ndofs))
end

function evaluate!(cache, k::MorleyChangeOfBasis, Jt, σ)
  nedges = length(k.tangents)
  nverts = nedges  # a triangle
  ndofs = nverts + nedges
  setsize!(cache, (ndofs, ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))

  detJ = det(Jt)
  G = inv(Jt ⋅ transpose(Jt))  # (JᵀJ)⁻¹, since Jt = Jᵀ

  for i in 1:nverts
    M[i, i] = 1.0
  end

  for e in 1:nedges
    t̂, n̂ = k.tangents[e], k.normals[e]
    Gn̂ = G ⋅ n̂
    A = detJ * (n̂ ⋅ Gn̂)
    B = detJ * (t̂ ⋅ Gn̂)
    σe = σ[e]
    va, vb = k.edge_vertices[e]

    # D = diag(1,1,1,σ...) scales the last columns of the block below.
    if k.transposed_inverse
      # WᵀD
      M[nverts+e, nverts+e] = A * σe
      M[va, nverts+e] = -B * σe
      M[vb, nverts+e] = B * σe
    else
      # W⁻¹D = [I 0; -A⁻¹B A⁻¹D]
      M[nverts+e, nverts+e] = σe / A
      M[nverts+e, va] = B / A
      M[nverts+e, vb] = -B / A
    end
  end

  return M
end

# DOF scaling: the default `h⁰` of the identity map is correct since the vertex
# values are invariant and the edge DoF `∫ₑ (∇u⋅n) ds` scales like `h ⋅ h⁻¹`.

