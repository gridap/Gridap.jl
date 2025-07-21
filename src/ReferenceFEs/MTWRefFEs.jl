
struct MardalTaiWinther <: ReferenceFEName end

const mtw = MardalTaiWinther()

Pushforward(::Type{MardalTaiWinther}) = ContraVariantPiolaMap()

"""
    struct MardalTaiWinther <: ReferenceFEName end
    MardalTaiWintherRefFE(::Type{et},p::Polytope,order::Integer) where et

Mardal-Tai-Winther reference finite element.

References:

- `A Robust Finite Element Method for Darcy-Stokes Flow`, Mardal, Tai and Winther (2002)
- `Transformations for Piola-mapped elements`, Aznaran, Farrell and Kirby (2022)

"""
function MardalTaiWintherRefFE(::Type{T},p::Polytope,order::Integer) where T
  D = num_dims(p)
  @assert is_simplex(p) "MardalTaiWinther Reference FE only defined simplices"
  @asset order == 3 "MardalTaiWinther Reference FE is by definition of order 3"
  # TODO: We should just not allow this to be an argument

  prebasis = MonomialBasis(Val(D),VectorValue{D,T},3,Polynomials._p_filter)
  eb = MonomialBasis(Val(1),T,0,Polynomials._p_filter)
  fb = MonomialBasis(Val(D-1),T,1,Polynomials._p_filter)

  function emom(φ,μ,ds) # Edge moment function: σ_K(φ,μ) = ∫((φ⋅t)*μ)dK
    t = get_edge_tangent(ds)
    φt = Broadcasting(Operation(⋅))(φ,t)
    Broadcasting(Operation(*))(φt,μ)
  end
  function fmom(φ,μ,ds) # Face moment function : σ_F(φ,μ) = ∫((φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end

  moments = [
    (get_dimrange(p,1),emom,eb),   # Edge moments
    (get_dimrange(p,D-1),fmom,fb), # Face moments
  ]

  return MomentBasedReferenceFE(MardalTaiWinther(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::MardalTaiWinther,::Type{T}, order) where T
  MardalTaiWintherRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{MardalTaiWinther},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a MardalTaiWinther reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{MardalTaiWinther}, conf::DivConformity)
  get_face_dofs(reffe)
end
