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
    HermiteRefFE(::Type{V}, p::Polytope; kwargs...)

The `kwargs` is [`poly_type`](@ref "`poly_type` keyword argument").
"""
function HermiteRefFE(
  ::Type{V}, p::Polytope{D}, order::Integer; poly_type=_mom_reffe_default_PT(p)) where {V,D}

  @check order == 3 "Hermite Reference FE only available for `order`=3, got order=$order"
  @check 1 ≤ D ≤ 3 && is_simplex(p) "Hermite Reference FE only available for simplices of dimension 1 to 3"

  PT = poly_type
  prebasis = FEEC_poly_basis(Val(D), V, order ,0,:P, PT) # PᵣΛᴰ⁻¹, r = order

  dofs = _hermite_dofs(V, p, prebasis)

  ndofs = length(prebasis)
  face_dofs = [Int[] for _ in 1:num_faces(p)]
  face_dofs[end] = 1:ndofs
  shapefuncs = compute_shapefuns(dofs, prebasis)
  conformity = C1Conformity()
  metadata = nothing

  return GenericRefFE{Bubble}(
    ndofs,
    p,
    prebasis,
    dofs,
    conformity,
    metadata,
    face_dofs,
    shapefuncs,
  )
end

"""
    _hermite_dofs(p::Polytope{D})

- at each vertex ``vᵢ`` of `p`: directional derivatives ``φ -> ∂ᵢφ(vᵢ)``
- at each barycenter ``bᵢ`` of facet ``fᵢ`` (``dim(fᵢ)>0``) of `p`: point-value ``φ -> φ(bᵢ)``
"""
function _hermite_dofs(::Type{V}, p::Polytope{D}, prebasis) where {V,D}

  vertices = get_vertex_coordinates(p)
  P = eltype(vertices)

  function partial_at_vertex(φ,μ,ds) # moment function: σ_K(φ,μ) = ∇φ⋅μ
    ∇φ = Broadcasting(∇)(φ) # this is what I do usually to eval gradient in RefFE
    #∇φ = Broadcasting(Operation(∇))(φ) # Does not work either
    Broadcasting(Operation(⋅))(∇φ,μ)
  end

  μ = MonomialBasis(Val(0),P,0) # cartesian vector basis
  moments = Tuple[ (get_dimrange(p,0), partial_at_vertex, μ), ]

  gradient_values = MomentBasedDofBasis(p, prebasis, moments)

  if D>1
    barycenters = Vector{P}(undef, num_facets(p))
    for (lfacet, facet_vertices) in enumerate(get_face_coordinates(p, D-1))
      barycenters[lfacet] = mean(facet_vertices)
    end
    bary_pointvalues = LagrangianDofBasis(V, barycenters)

    lazy_append(gradient_values, bary_pointvalues)
  else
    gradient_values
  end
end

function ReferenceFE(p::Polytope,::Hermite,::Type{V}, order; kwargs...) where V
  HermiteRefFE(V,p,order; kwargs...)
end

function Conformity(reffe::GenericRefFE{Hermite}, sym::Symbol)
  c1 = (:C1, :H2)
  if sym == :L2
    L2Conformity()
  elseif sym in c1
    C1Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Hermite reference FE.

    Possible values of conformity for this reference fe are $((:L2, c1...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{Hermite}, conf::C1Conformity)
  @warn "check get_face_own_dofs"
  get_face_dofs(reffe)
end
