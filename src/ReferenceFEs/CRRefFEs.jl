"""
    struct CR  <: ReferenceFEName
"""
struct CR <: ReferenceFEName end

"""
    const cr = CR()

Singleton of the [`CR`](@ref) reference FE name.
"""
const cr = CR()

"""
CRRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the P space of degree `order-1`.
"""
function CRRefFE(::Type{T},p::Polytope,order::Integer) where T
  D = num_dims(p)

  if is_simplex(p) && order == 1
    prebasis = MonomialBasis(Val(D),T,order,Polynomials._p_filter)
    fb = MonomialBasis(Val(D-1),T,0,Polynomials._p_filter)
  else
    @notimplemented "CR Reference FE only available for simplices and lowest order"
  end

  function fmom(φ,μ,ds) # Face moment function : σ_F(φ,μ) = 1/|F| ( ∫((φ)*μ)dF )
    D = num_dims(ds.cpoly)
    facet_measure = get_facet_measure(ds.cpoly, D-1)
    facet_measure_1 = Gridap.Fields.ConstantField(1 / facet_measure[ds.face])
    φμ = Broadcasting(Operation(⋅))(φ,μ)
    Broadcasting(Operation(*))(φμ,facet_measure_1)
  end

  moments = Tuple[
    (get_dimrange(p,D-1),fmom,fb), # Face moments
  ]

  return Gridap.ReferenceFEs.MomentBasedReferenceFE(CR(),p,prebasis,moments,L2Conformity())
end


function ReferenceFE(p::Polytope,::CR, order)
  CRRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::CR,::Type{T}, order) where T
  CRRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{CR},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a CR reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{CR}, conf::L2Conformity)
  get_face_dofs(reffe)
end
