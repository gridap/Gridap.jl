"""
    struct Argyris <: ReferenceFEName
"""
struct Argyris <: ReferenceFEName end

"""
    const argyris = Argyris()

Singleton of the [`Argyris`](@ref) reference FE name.
"""
const argyris = Argyris()

#Pushforward(::Type{Argyris}) = ContraVariantPiolaMap()

"""
    ArgyrisRefFE(::Type{V}, p::Polytope, order; vertices=nothing,
                 poly_type, mom_poly_type=poly_type)

The default conformity is `:H1`/`:Hgrad` because the global FE space is not C1,
although it is C1 at element vertices.

Available on simplices. If `V <: MultiValue`, a cartesian product of the scalar
Argyris FE is constructed.

The kwarg [`poly_type`](@ref "`poly_type` keyword argument") defaults to
`Bernstein`. `mom_poly_type` is used for the edge and face/cell moments for `order>5`.

If `vertices` are given, the DOF are defined on the physical triangle defined
by them.
"""
function ArgyrisRefFE(
  ::Type{V}, p::Polytope{D}, order; vertices=nothing,
  poly_type=_mom_reffe_default_PT(p), mom_poly_type=poly_type) where {V,D}

  @check D==2 && is_simplex(p) "Argyris Reference FE only available on triangles, got $p"
  @check 5 <= order "Argyris Reference FE only available for orders ≥ 5, got $order"
  #@notimplementedif order != 5

  PT = poly_type
  cart_prod = V <: MultiValue
  prebasis = FEEC_poly_basis(Val(D), V, order, 0,:P, PT; cart_prod, vertices) # PᵣΛ⁰(△₂)

  MT = mom_poly_type
  phys_p = isnothing(vertices) ? p : GeneralPolytope{D}(p, vertices)
  dofs, face_own_dofs = _argyris_dofs_and_faceowndofs(V, phys_p, prebasis, MT)

  ndofs = length(dofs)
  conformity = C1Conformity()
  metadata = nothing

  return GenericRefFE{Argyris}(
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
    _argyris_dofs_and_faceowndofs(V, p::Polytope{2}, prebasis, MT)

With ``r`` the order of the element (and `prebasis`), the DOF are:

Moments on φ
- at each vertex ``vᵢ`` of `p`: point-value ``φ -> φ(vᵢ)``
- at each edge ``eᵢ`` of `p`: moment ``φ -> ∫_f φ*μⱼ df`` for ``μⱼ ∈ ℙᵣ₋₆``
- at the face ``f`` of `p`: moment ``φ -> ∫_f φ*μⱼ df`` for ``μⱼ ∈ ℙᵣ₋₆``

Moments on ∇φ
- at each vertex ``vᵢ`` of `p`: directional derivatives ``φ -> ∂ⱼφ(vᵢ) = (∇φ)(vᵢ)ⱼ`` for ``1 ≤ j ≤ 2``
- at each edge ``eᵢ`` of `p`: normal derivative moment ``φ -> ∫_eᵢ (∇φ⋅n)*μⱼ deᵢ`` for ``μⱼ ∈ ℙᵣ₋₅``

Moments on ∇∇φ
- at each vertex ``vᵢ`` of `p`: directional derivatives ``φ -> ∂ₖ∂ⱼφ(vᵢ)`` for ``1 ≤ j ≤ k ≤ 2``
"""
function _argyris_dofs_and_faceowndofs(::Type{V}, p::Polytope{D}, prebasis, ::Type{MT}) where {V,D,MT<:Polynomial}
  P = eltype(get_vertex_coordinates(p))
  r = get_order(prebasis)
  cart_prod = V <: MultiValue


  # Μoments of φ

  # moment function: σ_K(φ,μ) = ∫_f φ⊙μ df (∫_f φ*μ df for scalar φ)
  function component_moment(φ,e,ds)
    Broadcasting(Operation(⊙))(φ,e)
  end

  vb = MonomialBasis(Val(0),V,0) # cartesian vector basis of V
  eb = r>5 ? FEEC_poly_basis(Val(1), V, r-6, 0, :P, MT; cart_prod) : nothing # Pᵣ₋₆Λ⁰(△₁) in each (indep) component of V
  cb = r>5 ? FEEC_poly_basis(Val(2), V, r-6, 0, :P, MT; cart_prod) : nothing # Pᵣ₋₆Λ⁰(△₂) in each (indep) component of V

  value_moments = Tuple[
    (get_dimrange(p,0), component_moment, vb), # φ point values at vertices
  ]
  r > 5 && push!(value_moments,
    (get_dimrange(p,1), component_moment, eb), # φ integral moments on edges
    (get_dimrange(p,2), component_moment, cb)  # φ integral moment on the face
  )
  value_dofs = MomentBasedDofBasis(p, prebasis, value_moments)


  # moment function: σ_K(φ,μ) =  ∇φ⊙μ (∇φ⋅μ for scalar φ) for gradient components
  #                       or    ∇∇φ⊙μ                     for hessian components
  function partials_at_vertex(∇φ,μ,ds)
    Broadcasting(Operation(⊙))(∇φ,μ)
  end

  # Μoments of ∇φ

  G = gradient_type(V,zero(P))
  grad_V_basis = [ p⊗v for p in component_basis(P) for v in component_basis(V)]
  grad_V_dual_basis = representatives_of_basis_dual(grad_V_basis)
  # Create a polynomial basis having constant basis polynomial equal to grad_V_basis
  e = MonomialBasis(Val(0),G,0)    # cartesian vector basis
  e_basis = evaluate(e, Point())   # actually same as component_basis(G)
  change = [ gd⊙eᵢ for eᵢ in e_basis, gd in grad_V_dual_basis ]
  vb_grad = linear_combination(change, e)
  @check grad_V_dual_basis == evaluate(vb_grad, Point())

  function normal_derivative_to_edge(∇φ,μ,ds) # σ_E(∇φ,μ) = ∫((n·∇φ)*μ)dE
    n = get_facet_normal(ds)
    ∂ₙφ = Broadcasting(Operation(⋅))(n,∇φ)
    Broadcasting(Operation(⊙))(∂ₙφ,μ)
  end
  eb_grad = FEEC_poly_basis(Val(1), V, r-5, 0, :P, MT; cart_prod) # Pᵣ₋₅Λ⁰(△₁)

  grad_moments = Tuple[
    (get_dimrange(p,0), partials_at_vertex, vb_grad),
    (get_dimrange(p,1), normal_derivative_to_edge, eb_grad),
  ]
  grad_dofs = MomentBasedDofBasis(p, prebasis, grad_moments, gradient)


  # Μoments of ∇∇φ

  T = eltype(V)
  Hs = SymTensorValue{D,T} # only consider the 3 independent 2nd order partial derivatives
  hess_V_basis = [ hi⊗v for hi in component_basis(Hs) for v in component_basis(V)]
  hess_V_dual_basis = representatives_of_basis_dual(hess_V_basis)

  H = gradient_type(G,zero(P))
  e = MonomialBasis(Val(0),H,0)    # cartesian vector basis
  e_basis = evaluate(e, Point())   # actually same as component_basis(H)
  change = [ hd⊙eᵢ for eᵢ in e_basis, hd in hess_V_dual_basis ]
  vb_hess = linear_combination(change, e)
  @check SArray.(hess_V_dual_basis) == SArray.(evaluate(vb_hess, Point()))

  hess_moments = Tuple[ (get_dimrange(p,0), partials_at_vertex, vb_hess), ]
  hess_dofs = MomentBasedDofBasis(p, prebasis, hess_moments, ∇∇)


  # Face own dofs, appending dof ids of the three bases

  value_fod = get_face_own_moments(value_dofs)
  grad_fod  = get_face_own_moments(grad_dofs)
  hess_fod  = get_face_own_moments(hess_dofs)

  face_own_dofs = deepcopy(value_fod)
  curr_dofs_number = length(value_dofs)
  for i in eachindex(face_own_dofs)
    append!(face_own_dofs[i], grad_fod[i] .+ curr_dofs_number)
  end
  curr_dofs_number += length(grad_dofs)
  for i in eachindex(face_own_dofs)
    append!(face_own_dofs[i], hess_fod[i] .+ curr_dofs_number)
  end

  dofs = vcat(value_dofs, grad_dofs, hess_dofs)
  dofs, face_own_dofs
end

function ReferenceFE(p::Polytope,::Argyris,::Type{V}, order; kwargs...) where V
  ArgyrisRefFE(V,p,order; kwargs...)
end

function Conformity(::GenericRefFE{Argyris}, sym::Symbol)
  c1 = (:H2, :C1)
  if sym == :L2
    L2Conformity()
  elseif sym in c1
    C1Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Argyris reference FE.

    Possible values of conformity for this reference fe are $((:L2, c1...)).
      """
  end
end

