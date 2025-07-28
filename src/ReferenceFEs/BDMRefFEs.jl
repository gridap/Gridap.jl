"""
    struct BDM <: ReferenceFEName
"""
struct BDM <: ReferenceFEName end

"""
    const bdm = BDM()

Singleton of the [`BDM`](@ref) reference FE name.
"""
const bdm = BDM()

Pushforward(::Type{BDM}) = ContraVariantPiolaMap()

"""
    BDMRefFE(::Type{T}, p::Polytope, order::Integer)

The `order` argument has the following meaning: the divergence of the  functions
in this basis is in the ℙ space of degree `order-1`. `T` is the type of scalar
components.
"""
function BDMRefFE(::Type{T},p::Polytope,order::Integer) where T
  @check order > 0 "BDM Reference FE only available for order > 0, got order=$order"
  D = num_dims(p)
  @check 2 ≤ D ≤ 3 && is_simplex(p) "BDM Reference FE only available for simplices of dimension 2 and 3"

  rotate_90 = (D==2)
  k = D-1
  PT = Bernstein # Could be a Kwargs, any basis works

  prebasis =     FEEC_poly_basis(Val(D),  T,order  ,k,:P, PT; rotate_90) # PᵣΛᴰ⁻¹
  fb =           FEEC_poly_basis(Val(D-1),T,order  ,0,:P⁻,Bernstein)            # Facet basis P⁻ᵨΛ⁰(△ᴰ⁻¹), ρ = r
  cb = order>1 ? FEEC_poly_basis(Val(D),  T,order-1,1,:P⁻,Bernstein) : nothing  # Cell basis  P⁻ᵨΛ¹(△ᴰ),   ρ = r-1

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
  if order > 1
    push!(moments,(get_dimrange(p,D),cmom,cb)) # Cell moments
  end

  return MomentBasedReferenceFE(BDM(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::BDM,::Type{T}, order) where T
  BDMRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{BDM},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a BDM reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{BDM}, conf::DivConformity)
  get_face_dofs(reffe)
end
