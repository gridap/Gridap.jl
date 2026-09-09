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

# Mapped by the plain pullback; this is also the default, stated for clarity.
Pushforward(::Type{Morley}) = IdentityPiolaMap()

# The Morley DoFs: the vertex point evaluations, then the edge normal-derivative
# moments. Both are `MomentBasedDofBasis`es -- the second built with `operator=∇`
# -- concatenated with `vcat`, since a single one carries a single operator for
# all of its moments.
function _morley_dof_basis(::Type{T}, p::Polytope{2}, prebasis) where T
  vb = MonomialBasis(Val(0), T, 0)  # the constant on a vertex
  eb = MonomialBasis(Val(1), T, 0)  # the constant on an edge

  vmom(φ, μ, ds) = Broadcasting(Operation(*))(φ, μ)      # σ_v(u,μ) = u(v) μ
  function emom(φ, μ, ds)                                # σ_e(∇u,μ) = ∫(∇u⋅n)μ
    n = _edge_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ, n)
    Broadcasting(Operation(*))(φn, μ)
  end

  vertex_dofs = MomentBasedDofBasis(p, prebasis, Tuple[(get_dimrange(p, 0), vmom, vb)])
  edge_dofs = MomentBasedDofBasis(p, prebasis, Tuple[(get_dimrange(p, 1), emom, eb)], ∇)
  return vcat(vertex_dofs, edge_dofs)
end

"""
    MorleyRefFE(::Type{T}, p::Polytope{2})

The Morley reference FE on the triangle `K`, with `T` the scalar type: the
quadratic nonconforming plate element of [Morley, Aero. Quart. 19 (1968) 149],
6 DoFs, with

    Morley(K) = P₂(K)

of dimension 6. The space needs no constraints; what is non-standard is the DoF
set, whose push-forward mixes the edge and vertex functionals — the element
therefore needs a cell-dependent change of basis, see `compute_cell_bases_changes`
in `FESpaces`.

## Prebasis

`P₂(K)`, of dimension 6, which is the element's space itself.

## Moments

At each vertex `v` the point value, and per edge `e` with unit tangent `t` and
normal `n = R t` (`R` the clockwise quarter turn) the mean normal derivative:

    ℓ^v(u)   = u(v)
    ℓ^e(u)   = ∫ₑ (∇u⋅n) ds

`3 + 3 = 6`. Writing the edge DoF as a *moment* rather than as the classical
midpoint value `∂u/∂n(m_e)` is what keeps the transformation exact. Under the
plain pullback `∇u⋅n` picks up a tangential derivative, which is not one of the
reference nodes, so the node completion of [Kirby, SMAI-JCM 4 (2018) 197, §3.3]
is called for; but in moment form that tangential part is
`∫ₑ (∇u⋅t) ds = u(v_b) - u(v_a)` by the fundamental theorem of calculus, for any
smooth `u` and not merely on `P₂`, so it is already a difference of two other
DoFs and the completion collapses. The two forms of the DoF agree up to the
factor `|e|` on `P₂` in any case, the normal derivative being linear along an
edge.

The conformity is `H1Conformity()`, which in Gridap fixes *which faces own and
share DoFs* — here all of them. It is not a claim that the assembled space is
C⁰: Morley is neither C⁰ nor C¹, the shape functions matching only at the
vertices and in the mean normal derivative across each edge.
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

function ReferenceFE(p::Polytope, ::Morley, ::Type{T}) where T
  MorleyRefFE(T, p)
end

function ReferenceFE(p::Polytope, ::Morley, ::Type{T}, order) where T
  @notimplementedif order != 2 """\n
  The Morley element exists for order 2 only, got $order.
  """
  MorleyRefFE(T, p)
end

# Identity for every admissible vertex permutation of every face: each face owns
# a single DoF, and reversing an edge changes only the sign of that DoF, which
# the change of basis carries.
function get_face_own_dofs_permutations(reffe::GenericRefFE{Morley}, conf::Conformity)
  _identity_dof_permutations(reffe, conf)
end
