"""
    struct Morley <: ReferenceFEName
"""
struct Morley <: ReferenceFEName end

"""
    const morley = Morley()

Singleton of the [`Morley`](@ref) reference FE name.
"""
const morley = Morley()

#Pushforward(::Type{Morley}) = ContraVariantPiolaMap()

"""
    MorleyRefFE(::Type{V}, p::Polytope, order; vertices=nothing,
                 poly_type, mom_poly_type=poly_type)

Available on simplices. If `V <: MultiValue`, a cartesian product of the scalar
Morley FE is constructed.

The default conformity is `:H1`/`:Hgrad` because the global FE space is not C1,
although shape function are normally continuous at edge midpoints.

The kwarg [`poly_type`](@ref "`poly_type` keyword argument") defaults to
`Bernstein`. `mom_poly_type` is used for the edge and face/cell moments for `order>5`.

If `vertices` are given, the DOF are defined on the physical triangle defined
by them.
"""
function MorleyRefFE(
  ::Type{V}, p::Polytope{D}, order; vertices=nothing,
  poly_type=_mom_reffe_default_PT(p)) where {V,D}

  @check D==2 && is_simplex(p) "Morley Reference FE only available on triangles, got $p"
  @check 2 == order "Morley Reference FE only available for orders 2, got $order"

  PT = poly_type
  cart_prod = V <: MultiValue
  prebasis = FEEC_poly_basis(Val(D), V, order, 0,:P, PT; cart_prod, vertices) # PᵣΛ⁰(△₂)

  phys_p = isnothing(vertices) ? p : GeneralPolytope{D}(p, vertices)
  dofs, face_own_dofs = _morley_dofs_and_faceowndofs(V, phys_p, prebasis)

  ndofs = length(dofs)
  conformity = GradConformity()
  metadata = nothing

  return GenericRefFE{Morley}(
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
    _morley_dofs_and_faceowndofs(V, p::Polytope{2}, prebasis)

With ``r`` the order of the element (and `prebasis`), the DOF are:

Moments on φ
- at each vertex ``vᵢ`` of `p`: point-value ``φ -> φ(vᵢ)``

Moments on ∇φ
- at each mid-point ``mᵢ`` of edge ``eᵢ`` of `p`: pointwise normal derivative ``φ -> (∇φ⋅n)|mᵢ``
"""
function _morley_dofs_and_faceowndofs(::Type{V}, p::Polytope{D}, prebasis) where {V,D}
  P = eltype(get_vertex_coordinates(p))
  cart_prod = V <: MultiValue
  # Μoments of φ

  # moment function: σ_K(φ,μ) = ∫_f φ⊙μ df (∫_f φ*μ df for scalar φ)
  function component_moment(φ,e,ds)
    Broadcasting(Operation(⊙))(φ,e)
  end

  vb = MonomialBasis(Val(0),V,0) # cartesian vector basis of V
  value_moments = [ (get_dimrange(p,0), component_moment, vb) ] # φ point values at vertices
  value_dofs = MomentBasedDofBasis(p, prebasis, value_moments)

  # Μoments of ∇φ

  # moment function: σ_K(φ,μ) =  ∇φ⊙μ (∇φ⋅μ for scalar φ) for gradient components

  G = gradient_type(V,zero(P))
  grad_V_basis = [ p⊗v for p in component_basis(P) for v in component_basis(V)]
  grad_V_dual_basis = representatives_of_basis_dual(grad_V_basis)
  # Create a polynomial basis having constant basis polynomial equal to grad_V_basis
  e = MonomialBasis(Val(1),G,0)    # cartesian vector basis
  e_basis = evaluate(e, Point(0.))   # actually same as component_basis(G)
  change = [ gd⊙eᵢ for eᵢ in e_basis, gd in grad_V_dual_basis ]
  mi_grad = linear_combination(change, e)
  @check grad_V_dual_basis == evaluate(mi_grad, Point(0.))

  mi_grad = FEEC_poly_basis(Val(1), V, 0, 0, :P; cart_prod) # Pᵣ₋₅Λ⁰(△₁)

  function normal_derivative_to_edge(∇φ,μ,ds) # σ_E(∇φ,μ) = μ*(n·∇φ)|mᵢ
    n = get_facet_normal(ds)
    ∂ₙφ = Broadcasting(Operation(⋅))(n,∇φ)

    # only necessary as long as PointValueDofBasis does not cancel ref-face volume
    face_measure = ReferenceFEs._get_dfaces_measure(ds.cpoly, 1)
    face_measure_1 = Gridap.Fields.ConstantField(1 / getindex(face_measure,ds.face))
    detJ_∂ₙφ  = Broadcasting(Operation(*))(∂ₙφ,face_measure_1)

    Broadcasting(Operation(⊙))(detJ_∂ₙφ,μ)
  end

  #μ = MonomialBasis{1}(Float64,0)
  mid_SEGMENT_point = mean(get_vertex_coordinates(SEGMENT))
  grad_moments = Tuple[ (get_dimrange(p,1), normal_derivative_to_edge, mi_grad, [mid_SEGMENT_point]) ]
  grad_dofs = PointValueDofBasis(p, prebasis, grad_moments, gradient)

  # Face own dofs, appending dof ids of the three bases

  value_fod = get_face_own_moments(value_dofs)
  grad_fod  = get_face_own_moments(grad_dofs)

  face_own_dofs = deepcopy(value_fod)
  curr_dofs_number = length(value_dofs)
  for i in eachindex(face_own_dofs)
    append!(face_own_dofs[i], grad_fod[i] .+ curr_dofs_number)
  end

  dofs = vcat(value_dofs, grad_dofs)
  dofs, face_own_dofs
end

function ReferenceFE(p::Polytope,::Morley,::Type{V}, order; kwargs...) where V
  MorleyRefFE(V,p,order; kwargs...)
end

function Conformity(::GenericRefFE{Morley}, sym::Symbol)
  h1 = (:H1, :C0)
  if sym == :L2
    L2Conformity()
  elseif sym in h1
    GradConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Morley reference FE.

    Possible values of conformity for this reference fe are $((:L2, h1...)).
      """
  end
end

