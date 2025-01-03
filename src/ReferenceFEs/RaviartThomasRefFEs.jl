
struct DivConformity <: Conformity end
abstract type DivConforming <: ReferenceFEName end

# RaviartThomas

struct RaviartThomas <: DivConforming end
const raviart_thomas = RaviartThomas()

"""
    RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function RaviartThomasRefFE(
  ::Type{T},p::Polytope{D},order::Integer
) where {T,D}

  if is_n_cube(p)
    prebasis = QCurlGradBasis(Legendre,Val(D),T,order)         # Prebasis
    cb = QGradBasis(Legendre,Val(D),T,order-1)                 # Cell basis
    fb = LegendreBasis(Val(D-1),T,order,Polynomials._q_filter) # Face basis
  elseif is_simplex(p)
    prebasis = PCurlGradBasis(Monomial,Val(D),T,order)                        # Prebasis
    cb = LegendreBasis(Val(D),VectorValue{D,T},order-1,Polynomials._p_filter) # Cell basis
    fb = LegendreBasis(Val(D-1),T,order,Polynomials._p_filter)                # Face basis
  else
    @notimplemented "Raviart-Thomas Reference FE only available for cubes and simplices"
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

function ReferenceFE(p::Polytope,::RaviartThomas,order;kwargs...)
  RaviartThomasRefFE(Float64,p,order;kwargs...)
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

# TODO: Please remove me
function legendreBasis(::Type{T},p::Polytope,orders) where T
  compute_legendre_basis(T,p,orders)
end
function legendreBasis(::Type{T},p::Polytope{D},order::Int) where {D,T}
  orders = tfill(order,Val{D}())
  legendreBasis(T,p,orders)
end
function compute_legendre_basis(::Type{T},p::ExtrusionPolytope{D},orders) where {D,T}
  extrusion = Tuple(p.extrusion)
  terms = _monomial_terms(extrusion,orders)
  LegendreBasis(Val(D),T,orders,terms)
end

# ContraVariantPiolaMap

struct ContraVariantPiolaMap <: Map end

function evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Fields.BroadcastOpFieldArray{ContraVariantPiolaMap}
)
  v, Jt, sign_flip = a.args
  ∇v = Broadcasting(∇)(v)
  k = ContraVariantPiolaMap()
  Broadcasting(Operation(k))(∇v,Jt,sign_flip)
end

function lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::LazyArray{<:Fill{Broadcasting{Operation{ContraVariantPiolaMap}}}}
)
  v, Jt, sign_flip = a.args
  ∇v = lazy_map(Broadcasting(∇),v)
  k = ContraVariantPiolaMap()
  lazy_map(Broadcasting(Operation(k)),∇v,Jt,sign_flip)
end

function lazy_map(
  k::ContraVariantPiolaMap,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field},
  sign_flip::AbstractArray{<:AbstractArray{<:Field}}
)
  cell_Jt = lazy_map(∇,cell_map)
  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_Jt,sign_flip)
end

function evaluate!(
  cache,::ContraVariantPiolaMap,
  v::Number,
  Jt::Number,
  sign_flip::Bool
)
  idetJ = 1/meas(Jt)
  ((-1)^sign_flip*v)⋅(idetJ*Jt)
end

function evaluate!(
  cache,
  k::ContraVariantPiolaMap,
  v::AbstractVector{<:Field},
  phi::Field,
  sign_flip::AbstractVector{<:Field}
)
  Jt = ∇(phi)
  Broadcasting(Operation(k))(v,Jt,sign_flip)
end
