
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
    RaviartThomasRefFE(::Type{T}, p::Polytope, order::Integer)

The `order` argument has the following meaning: the divergence of the functions
in this basis is in the Q space of degree `order`. `T` is the type of scalar components.
"""
function RaviartThomasRefFE(
  ::Type{T},p::Polytope{D},order::Integer
) where {T,D}

  rotate_90 = D==2
  k = D-1

  if is_n_cube(p)
    PT = Legendre # Could be a Kwargs, any hierarchical basis works
    @check PT ≠ Bernstein # broken for Bernstein in prebasis or cb, might be an issue of ordering of the basis polynomials
    prebasis =     FEEC_poly_basis(Val(D),  T,order+1,k,:Q⁻,PT; rotate_90) # Q⁻ᵣΛᵏ(□ᴰ), r = order+1
    fb =           FEEC_poly_basis(Val(D-1),T,order  ,0,:Q⁻,PT)            # Facet basis Q⁻ᵨΛ⁰(□ᴰ⁻¹), ρ = r-1
    cb = order>0 ? FEEC_poly_basis(Val(D),  T,order  ,1,:Q⁻,PT) : nothing  # Cell basis  Q⁻ᵨΛ¹(□ᴰ),   ρ = r-1
  elseif is_simplex(p)
    PT = Bernstein # Could be a Kwargs, any basis works
    prebasis =     FEEC_poly_basis(Val(D),  T,order+1,k,:P⁻,PT; rotate_90) # P⁻ᵣΛᵏ(△ᴰ), r = order+1
    fb =           FEEC_poly_basis(Val(D-1),T,order  ,0,:P ,PT)            # Facet basis PᵨΛ⁰(△ᴰ⁻¹), ρ = r-1
    cb = order>0 ? FEEC_poly_basis(Val(D),  T,order-1,1,:P ,PT) : nothing  # Cell basis  PᵨΛ¹(△ᴰ),   ρ = r-2
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

  return MomentBasedReferenceFE(RaviartThomas(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::RaviartThomas,::Type{T},order;kwargs...) where T
  RaviartThomasRefFE(T,p,order;kwargs...)
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
