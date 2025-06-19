"""
    struct BDM <: PushforwardRefFE <: ReferenceFEName
"""
struct BDM <: PushforwardRefFE end

"""
    const bdm = BDM()

Singleton of the [`BDM`](@ref) reference FE name.
"""
const bdm = BDM()

Pushforward(::Type{<:BDM}) = ContraVariantPiolaMap()

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

  prebasis = MonomialBasis(Val(D),VectorValue{D,T},order,Polynomials._p_filter)
  fb = MonomialBasis(Val(D-1),T,order,Polynomials._p_filter)
  cb = PGradBasis(Monomial,Val(D),T,order-2)

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
