
# RaviartThomas

"""
    struct RaviartThomas <: ReferenceFEName
"""
struct RaviartThomas <: ReferenceFEName end

"""
    const raviart_thomas = RaviartThomas()

Singleton of the [`RaviartThomas`](@ref) reference FE name.
"""
const raviart_thomas = RaviartThomas()

Pushforward(::Type{RaviartThomas}) = ContraVariantPiolaMap()

"""
    RaviartThomasRefFE(::Type{T}, p::Polytope, order::Integer; kwargs...)

The `order` argument has the following meaning: the divergence of the functions
in this basis is in the Q space of degree `order`. `T` is the type of scalar components.

The `kwargs` are [`change_dof`](@ref "`change_dof` keyword argument"),
[`poly_type`](@ref "`poly_type` keyword argument") and
[`mom_poly_type`](@ref "`mom_poly_type` keyword argument").
"""
function RaviartThomasRefFE(
  ::Type{T},p::Polytope{D},order::Integer;
  change_dof=true, poly_type=_mom_reffe_default_PT(p), mom_poly_type=poly_type) where {T,D}

  PT, MPT = poly_type, mom_poly_type
  rotate_90 = D==2
  k = D-1

  if is_n_cube(p)
    prebasis =     FEEC_poly_basis(Val(D),  T,order+1,k,:Q⁻,PT; rotate_90) # Q⁻ᵣΛᵏ(□ᴰ), r = order+1
    fb =           FEEC_poly_basis(Val(D-1),T,order  ,0,:Q⁻,MPT)            # Facet basis Q⁻ᵨΛ⁰(□ᴰ⁻¹), ρ = r-1
    cb = order>0 ? FEEC_poly_basis(Val(D),  T,order  ,1,:Q⁻,MPT) : nothing  # Cell basis  Q⁻ᵨΛ¹(□ᴰ),   ρ = r-1
  elseif is_simplex(p)
    prebasis =     FEEC_poly_basis(Val(D),  T,order+1,k,:P⁻,PT; rotate_90) # P⁻ᵣΛᵏ(△ᴰ), r = order+1
    fb =           FEEC_poly_basis(Val(D-1),T,order  ,0,:P ,MPT)            # Facet basis PᵨΛ⁰(△ᴰ⁻¹), ρ = r-1
    cb = order>0 ? FEEC_poly_basis(Val(D),  T,order-1,1,:P ,MPT) : nothing  # Cell basis  PᵨΛ¹(△ᴰ),   ρ = r-2
  else
    @notimplemented "Raviart-Thomas Reference FE only available for n-cubes and simplices"
  end

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function : σ_F(φ,μ) = ∫((φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end

  moments = Tuple[
    (get_dimrange(p,D-1),fmom,fb), # Face moments
  ]
  if (order > 0)
    push!(moments,(get_dimrange(p,D),cmom,cb)) # Cell moments
  end

  conf = DivConformity()
  change_dof = _validate_change_dof(change_dof, prebasis, p, conf)
  return MomentBasedReferenceFE(raviart_thomas,p,prebasis,moments,conf; change_dof)
end

function ReferenceFE(p::Polytope,::RaviartThomas,::Type{T},order; kwargs...) where T
  RaviartThomasRefFE(T,p,order; kwargs...)
end

function Conformity(reffe::GenericRefFE{RaviartThomas},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Raviart Thomas reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
    """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{RaviartThomas}, conf::DivConformity)
  get_face_dofs(reffe)
end
