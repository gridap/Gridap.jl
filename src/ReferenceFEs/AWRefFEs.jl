struct ArnoldWinther <: ReferenceFEName end

const arnoldwinther = ArnoldWinther()

"""
ArnoldWintherRefFE(::Type{T},p::Polytope,order::Integer) where T

Arnold-Winther reference finite element.

References:

- `Mixed Finite Elements for Elasticity`, Arnold and Winther (2002)
- `Nonconforming Mixed Finite Elements for Elasticity`, Arnold and Winther (2003)
- `Transformations for Piola-mapped elements`, Aznaran, Farrell and Kirby (2022)

"""
function ArnoldWintherRefFE(::Type{T},p::Polytope,order::Integer) where T
  @assert p == TRI "ArnoldWinther Reference FE only defined for TRIangles"
  conforming = true # TODO: Make this an argument
  
  VT = SymTensorValue{2,T}
  prebasis = MonomialBasis{2}(VT,3,Polynomials._p_filter)
  fb = MonomialBasis{D-1}(T,0,Polynomials._p_filter)
  cb = map(constant_field,component_basis(VT))

  function cmom(φ,μ,ds) # Cell and Node moment function: σ_K(φ,μ) = ∫(φ:μ)dK
    Broadcasting(Operation(⊙))(φ,μ)
  end
  function fmom_n(φ,μ,ds) # Face moment function (normal) : σ_F(φ,μ) = ∫((n·φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    nφn = Broadcasting(Operation(⋅))(n,φn)
    Broadcasting(Operation(*))(nφn,μ)
  end
  function fmom_t(φ,μ,ds) # Face moment function (tangent) : σ_F(φ,μ) = ∫((n·φ·t)*μ)dF
    n = get_facet_normal(ds)
    t = get_edge_tangent(ds)
    φn = Broadcasting(Operation(⋅))(φ,t)
    nφn = Broadcasting(Operation(⋅))(n,φn)
    Broadcasting(Operation(*))(nφn,μ)
  end

  moments = Tuple[
    (get_dimrange(p,1),fmom_n,fb), # Face moments (normal-normal)
    (get_dimrange(p,1),fmom_t,fb), # Face moments (normal-tangent)
    (get_dimrange(p,2),cmom,cb)    # Cell moments
  ]

  if conforming
    node_moments = Tuple[(get_dimrange(p,0),cmom,cb)]  # Node moments
    moments = vcat(node_moments,moments)
  end

  return MomentBasedReferenceFE(ArnoldWinther(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::ArnoldWinther, order)
  BDMRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::ArnoldWinther,::Type{T}, order) where T
  BDMRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{ArnoldWinther},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a ArnoldWinther reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{ArnoldWinther}, conf::DivConformity)
  get_face_dofs(reffe)
end
