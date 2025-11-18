"""
    struct Hermite <: ReferenceFEName
"""
struct Hermite <: ReferenceFEName end

"""
    const hermite = Hermite()

Singleton of the [`Hermite`](@ref) reference FE name.
"""
const hermite = Hermite()

#Pushforward(::Type{Hermite}) = ContraVariantPiolaMap()

"""
    HermiteRefFE(::Type{V}, p::Polytope; poly_type)

The default conformity is `:H1`/`:Hgrad` because the global FE space is not C1,
although it is C1 at element vertices.

Available on simplices. If `V <: MultiValue`, a cartesian product of the scalar
Hermite FE is constructed.

The kwarg [`poly_type`](@ref "`poly_type` keyword argument") defaults to
`Bernstein`.
"""
function HermiteRefFE(
  ::Type{V}, p::Polytope{D}; poly_type=_mom_reffe_default_PT(p)) where {V,D}

  @check 1 ≤ D && is_simplex(p) "Hermite Reference FE only available on simplices of dimension ≥ 1"

  PT = poly_type
  cart_prod = V <: MultiValue
  prebasis = FEEC_poly_basis(Val(D), V, 3 ,0,:P, PT; cart_prod) # P₃Λᴰ⁻¹

  dofs, face_own_dofs = _hermite_dofs_and_faceowndofs(V, p, prebasis)

  ndofs = length(dofs)
  conformity = GradConformity()
  metadata = nothing

  return GenericRefFE{Hermite}(
    ndofs,
    p,
    prebasis,
    dofs,
    conformity,
    metadata,
    face_own_dofs,
  )
end

"""
    _hermite_dofs_and_faceowndofs(p::Polytope{D})

- at each vertex ``vᵢ`` of `p`: directional derivatives ``φ -> ∂ᵢφ(vᵢ)``
- at each vertex ``vᵢ`` of `p`: point-value ``φ -> φ(vᵢ)``
- at each barycenter ``bᵢ`` of facet ``fᵢ`` (``dim(fᵢ)>0``) of `p`: point-value ``φ -> φ(bᵢ)``
"""
function _hermite_dofs_and_faceowndofs(::Type{V}, p::Polytope{D}, prebasis) where {V,D}
  P = eltype(get_vertex_coordinates(p))

  G = gradient_type(V,zero(P))
  grad_V_basis = [ p⊗v for p in component_basis(P) for v in component_basis(V)]
  grad_V_dual_basis = representatives_of_basis_dual(grad_V_basis)

  # Create a polynomial basis having constant basis polynomial equal to grad_V_basis
  e = MonomialBasis(Val(0),G,0)    # cartesian vector basis
  e_basis = evaluate(e, Point())   # actually same as component_basis(G)
  change = [ gd⊙eᵢ for eᵢ in e_basis, gd in grad_V_dual_basis ]
  μ = linear_combination(change, e)
  @check grad_V_dual_basis == evaluate(μ, Point())

  # moment function: σ_K(φ,μ) = ∇φ⊙μ (∇φ⋅μ for scalar φ)
  function partials_at_vertex(∇φ,e,ds)
    Broadcasting(Operation(⊙))(∇φ,e)
  end
  moments = Tuple[ (get_dimrange(p,0), partials_at_vertex, μ), ]
  partials_pointvalues = MomentBasedDofBasis(p, prebasis, moments, gradient)

  nodes = P[]
  node_nb = 1
  face_own_nodes = [ Int[] for _ in 1:num_faces(p) ]
  for (vertex, vertex_coordinates) in zip(get_dimrange(p,0), get_face_coordinates(p, 0))
    push!(nodes, first(vertex_coordinates))
    push!(face_own_nodes[vertex], node_nb)
    node_nb += 1
  end
  if D≥2
    for (face, face_coordinates) in zip(get_dimrange(p, 2), get_face_coordinates(p, 2))
      push!(nodes, mean(face_coordinates))
      push!(face_own_nodes[face], node_nb)
      node_nb += 1
    end
  end
  pointvalues = LagrangianDofBasis(V, nodes)

  face_own_dofs = _generate_face_own_dofs(face_own_nodes, pointvalues.node_and_comp_to_dof)
  pp_face_own_dofs = get_face_own_moments(partials_pointvalues)
  # appending these two face_own_dofs
  nb_nodal_dofs = length(pointvalues)
  for i in eachindex(face_own_dofs)
    append!(face_own_dofs[i], pp_face_own_dofs[i] .+ nb_nodal_dofs)
  end

  dofs = vcat(pointvalues, partials_pointvalues)
  dofs, face_own_dofs
end

function ReferenceFE(p::Polytope,::Hermite,::Type{V}, order; kwargs...) where V
  @check order == 3 "Hermite Reference FE only available for `order`=3, got order=$order"
  HermiteRefFE(V,p; kwargs...)
end

function Conformity(::GenericRefFE{Hermite}, sym::Symbol)
  hgrad = (:H1, :C0)
  if sym == :L2
    L2Conformity()
  elseif sym in hgrad
    GradConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Hermite reference FE.

    Possible values of conformity for this reference fe are $((:L2, hgrad...)).
      """
  end
end

