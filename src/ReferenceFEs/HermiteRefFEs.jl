"""
    struct Hermite <: ReferenceFEName end

Reference FE name for the cubic Hermite simplex, a C⁰ element with C¹
continuity at the vertices. See [`HermiteRefFE`](@ref). Its singleton instance
is [`hermite`](@ref).
"""
struct Hermite <: ReferenceFEName end

"""
    const hermite = Hermite()

Singleton of the [`Hermite`](@ref) reference FE name.
"""
const hermite = Hermite()

Pushforward(::Type{Hermite}) = IdentityPiolaMap()

"""
    HermiteRefFE(::Type{V}, K::Polytope{D}; vertices=nothing, poly_type)

The cubic Hermite reference FE on the simplex `K` of dimension `D ≥ 1`, with
`V` the value type [Ciarlet & Raviart, Arch. Rational Mech. Anal. 46 (1972) 177].
Its conformity is `:H1`: the global space is C⁰, and C¹ at the vertices only.

The implementation follows [Kirby, SMAI-JCM 4 (2018) 197, §3.2].

If `V <: MultiValue`, the Cartesian product of `num_indep_components(V)` copies
of the scalar element is built, with the value type stacked into the shape
functions and the DoFs.

The kwarg [`poly_type`](@ref "`poly_type` keyword argument") defaults to
`Bernstein`.

If `vertices` are given, the DoFs and the shape functions are those of the
physical simplex with these vertices instead of the reference one, so the
element is defined directly on that cell.

# Extended help

## Prebasis

The prebasis is taken as `P₃(K;V)`, of dimension `(D+1)(D+2)(D+3)/6 ⋅ num_indep_components(V)`.

## DoFs

At each vertex `v` the value and the `D` Cartesian gradient components, and at
the barycenter `b` of each 2-face (the cell itself for `D = 2`) the value:

    ℓ^{v,0}(u) = u(v)
    ℓ^{v,i}(u) = ∇u(v) ⋅ eᵢ,    i = 1 … D
    ℓ^{b}(u)   = u(b)

for `V` scalar, and the same functionals on every component otherwise.
"""
function HermiteRefFE(::Type{V}, p::Polytope{D};
                      vertices=nothing, poly_type=_mom_reffe_default_PT(p)) where {V,D}
  @notimplementedif !(1 ≤ D && is_simplex(p)) """\n
  The Hermite element is only defined on simplices of dimension ≥ 1, got $p.
  """
  cart_prod = V <: MultiValue
  prebasis = FEEC_poly_basis(Val(D), V, 3, 0, :P, poly_type; cart_prod, vertices)  # P₃Λ⁰(△_D)
  dofs, face_own_dofs = _hermite_dofs_and_face_own_dofs(V, p, prebasis, vertices)

  ndofs = length(dofs)
  @check ndofs == length(prebasis)

  GenericRefFE{Hermite}(ndofs, p, prebasis, dofs, H1Conformity(), nothing, face_own_dofs)
end

function ReferenceFE(p::Polytope, ::Hermite, ::Type{V}; kwargs...) where V
  HermiteRefFE(V, p; kwargs...)
end

function ReferenceFE(p::Polytope, ::Hermite, ::Type{V}, order; kwargs...) where V
  @notimplementedif order != 3 """\n
  The Hermite element exists for order 3 only, got $order.
  """
  HermiteRefFE(V, p; kwargs...)
end

function _hermite_dofs_and_face_own_dofs(::Type{V}, p::Polytope{D}, prebasis, vertices) where {V,D}
  P = eltype(get_vertex_coordinates(p))

  nodes = P[]
  face_own_nodes = [Int[] for _ in 1:num_faces(p)]
  for (vertex, vertex_coordinates) in zip(get_dimrange(p, 0), get_face_coordinates(p, 0))
    push!(nodes, first(vertex_coordinates))
    push!(face_own_nodes[vertex], length(nodes))
  end
  if D ≥ 2
    for (face, face_coordinates) in zip(get_dimrange(p, 2), get_face_coordinates(p, 2))
      push!(nodes, mean(face_coordinates))
      push!(face_own_nodes[face], length(nodes))
    end
  end
  # The nodal values are a lagrangian nodal basis because we need dofs at the
  # barycenters of the 2-faces.
  values = LagrangianDofBasis(V, _hermite_map_nodes(nodes, vertices))

  grad_V_basis = [eᵢ ⊗ vⱼ for eᵢ in component_basis(P) for vⱼ in component_basis(V)]
  μ = map(constant_field, representatives_of_basis_dual(grad_V_basis))
  Dmom(∇φ, μ, ds) = Broadcasting(Operation(⊙))(∇φ, μ)  # ∇u(v) ⊙ μ
  grads = MomentBasedDofBasis(p, prebasis, Tuple[(get_dimrange(p, 0), Dmom, μ)], ∇)
  # move the moment quadrature points on the triangle with `vertice` vertex coordinates, if not nothing
  grads = MomentBasedDofBasis(
    _hermite_map_nodes(get_nodes(grads), vertices), get_face_moments(grads),
    get_face_nodes_dofs(grads), get_face_own_moments(grads), ∇
  )

  face_own_dofs = _generate_face_own_dofs(face_own_nodes, values.node_and_comp_to_dof)
  nvalues = length(values)
  for (own, own_grads) in zip(face_own_dofs, get_face_own_moments(grads))
    append!(own, own_grads .+ nvalues)
  end

  return vcat(values, grads), face_own_dofs
end

_hermite_map_nodes(nodes, ::Nothing) = nodes

# the affine map of the reference simplex onto `vertices`, vertex to vertex
function _hermite_map_nodes(nodes, vertices::AbstractVector{<:Point{D}}) where D
  v₁ = first(vertices)
  map(nodes) do ξ
    v₁ + sum(i -> (vertices[i+1] - v₁) * ξ[i], 1:D)
  end
end

_hermite_map_nodes(nodes, vertices::Tuple) = _hermite_map_nodes(nodes, collect(vertices))

# every DoF is a point evaluation in the global Cartesian frame at a point
# that any relabeling of the face's vertices fixes -- a vertex, or a barycenter
function get_face_own_dofs_permutations(reffe::GenericRefFE{Hermite}, conf::Conformity)
  _identity_dof_permutations(reffe, conf)
end

################################################################################
# Change of basis
#
# Hermite is mapped by the plain pullback u = û∘F⁻¹. The values are preserved:
# F is affine on a simplex, so it sends the reference vertices and barycenters to
# the physical ones. The vertex gradients are not: with K = J⁻ᵀ,
#
#   ∇u(v) = K ∇û(v̂),
#
# so the push-forward of the Cartesian derivatives at a vertex is a combination
# of the reference ones at the same vertex -- Kirby's affine-interpolation
# equivalence, (3.15) -- and
#
#   W = blockdiag(I, Kg, …, Kg),   W⁻¹ = blockdiag(I, Kg⁻¹, …, Kg⁻¹),
#
# one block Kg = K ⊗ I per vertex over its gradient DoFs, which are ordered
# direction by direction with the value component running fastest, so K acts on
# the direction index alone. Kg⁻¹ is the same construction on K⁻¹ = Jᵀ, so no
# matrix is inverted numerically for P = W⁻¹; P⁻ᵀ = Wᵀ needs K itself.
#
# There is no orientation sign flips because all DoFs are point valued.
#
# Nothing above depends on the source simplex being the reference one: `Jt` is
# the (transposed) Jacobian of the affine map from the simplex the element's
# DoFs are stated on -- `vertices`, if given -- onto the target simplex, vertex
# to vertex in order.
struct HermiteChangeOfBasis <: Map
  vertex_grad_dofs::Vector{Vector{Int}}  # per vertex, direction-major
  ndofs::Int
  transposed_inverse::Bool
end

function HermiteChangeOfBasis(reffe::GenericRefFE{Hermite}, transposed_inverse::Bool)
  HermiteChangeOfBasis(_hermite_vertex_grad_dofs(reffe), num_dofs(reffe), transposed_inverse)
end

# per vertex, the gradient DoFs, direction-major
function _hermite_vertex_grad_dofs(reffe::GenericRefFE{Hermite})
  p = get_polytope(reffe)
  D = num_dims(p)
  own = get_face_own_dofs(reffe)
  # a vertex owns its value(s) first, then D gradient components per value
  map(get_dimrange(p, 0)) do v
    nc = length(own[v]) ÷ (D + 1)
    own[v][nc+1:end]
  end
end

function return_cache(k::HermiteChangeOfBasis, Jt)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::HermiteChangeOfBasis, Jt)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0     # the values are preserved
  end

  # W carries K = J⁻ᵀ on the gradients, so W⁻¹ carries K⁻¹ = Jᵀ = Jt while Wᵀ
  # carries Kᵀ.
  Kg = k.transposed_inverse ? transpose(inv(Jt)) : Jt
  D = size(Kg, 1)

  for dofs in k.vertex_grad_dofs
    nc = length(dofs) ÷ D
    for c in 1:nc, i in 1:D, j in 1:D
      M[dofs[(i-1)*nc+c], dofs[(j-1)*nc+c]] = Kg[i, j]
    end
  end

  return M
end

################################################################################
# DOF scaling
#
# The point values are invariant, `h⁰`. A Cartesian derivative at a vertex
# carries one `K = J⁻ᵀ`, so it scales like `h⁻¹`.
function get_dofscale_setter_function(reffe::GenericRefFE{Hermite}, ::IdentityPiolaMap)
  exponent = zeros(Int, num_dofs(reffe))
  for dofs in _hermite_vertex_grad_dofs(reffe)
    exponent[dofs] .= -1
  end
  _dofscale_setter_from_exponents(exponent)
end

