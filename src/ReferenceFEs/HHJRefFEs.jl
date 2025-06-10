
struct HellanHerrmannJhonson <: ReferenceFEName end

const hhj = HellanHerrmannJhonson()

Pushforward(::Type{HellanHerrmannJhonson}) = DoubleContraVariantPiolaMap()

"""
    struct HellanHerrmannJhonson <: ReferenceFEName end
    HellanHerrmannJhonsonRefFE(::Type{T},p::Polytope,order::Integer) where T

Hellan-Herrmann-Jhonson reference finite element.

References:

- `The Hellan-Herrmann-Johnson method with curved elements`, Arnold and Walker (2020)

"""
function HellanHerrmannJhonsonRefFE(::Type{T},p::Polytope,order::Integer) where T
  @assert p == TRI "HellanHerrmannJhonson Reference FE only defined for TRIangles"

  VT = SymTensorValue{2,T}
  prebasis = MonomialBasis{2}(VT,order,Polynomials._p_filter)
  fb = MonomialBasis{D-1}(T,order,Polynomials._p_filter)
  cb = MonomialBasis{2}(VT,order-1,Polynomials._p_filter)

  function cmom(φ,μ,ds) # Cell and Node moment function: σ_K(φ,μ) = ∫(φ:μ)dK
    Broadcasting(Operation(⊙))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function (normal) : σ_F(φ,μ) = ∫((n·φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    nφn = Broadcasting(Operation(⋅))(n,φn)
    Broadcasting(Operation(*))(nφn,μ)
  end

  moments = Tuple[
    (get_dimrange(p,1),fmom,fb), # Face moments
    (get_dimrange(p,2),cmom,cb)  # Cell moments
  ]

  return MomentBasedReferenceFE(HellanHerrmannJhonson(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::HellanHerrmannJhonson, order)
  HellanHerrmannJhonsonRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::HellanHerrmannJhonson,::Type{T}, order) where T
  HellanHerrmannJhonsonRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{HellanHerrmannJhonson},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a HellanHerrmannJhonson reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{HellanHerrmannJhonson}, conf::DivConformity)
  get_face_dofs(reffe)
end
