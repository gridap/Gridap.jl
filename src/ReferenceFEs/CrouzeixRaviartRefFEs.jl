"""
    struct CrouzeixRaviart  <: ReferenceFEName
"""
struct CrouzeixRaviart <: ReferenceFEName end

"""
    const couzeix_raviart = CrouzeixRaviart()

Singleton of the [`CrouzeixRaviart`](@ref) reference FE name.
"""
const crouzeix_raviart = CrouzeixRaviart()

"""
    CrouzeixRaviartRefFE(::Type{T}, p::Polytope, order::Integer)

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the ℙ space of degree `order-1`. `T` is the type of scalar components.
"""
function CrouzeixRaviartRefFE(::Type{T},p::Polytope,order::Integer) where T
  D = num_dims(p)

  if is_simplex(p) && order == 1
    prebasis = MonomialBasis(Val(D),T,order,Polynomials._p_filter)
    fb = MonomialBasis(Val(D-1),T,0,Polynomials._p_filter)
  else
    @notimplemented "Crouzet-Raviart Reference FE only available for simplices and lowest order"
  end

  function fmom(φ,μ,ds) # Face moment function : σ_F(φ,μ) = 1/|F| ( ∫((φ)*μ)dF )
    D = num_dims(ds.cpoly)
    facet_measure = _get_dfaces_measure(ds.cpoly, D-1)
    facet_measure_1 = Gridap.Fields.ConstantField(1 / facet_measure[ds.face])
    φμ = Broadcasting(Operation(⋅))(φ,μ)
    Broadcasting(Operation(*))(φμ,facet_measure_1)
  end

  moments = Tuple[
    (get_dimrange(p,D-1),fmom,fb), # Face moments
  ]

  return Gridap.ReferenceFEs.MomentBasedReferenceFE(CrouzeixRaviart(),p,prebasis,moments,L2Conformity())
end

function ReferenceFE(p::Polytope,::CrouzeixRaviart,::Type{T}, order) where T
  CrouzeixRaviartRefFE(T,p,order)
end


function Conformity(reffe::GenericRefFE{CrouzeixRaviart},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a CrouzeixRaviart reference FE.

    Possible values of conformity for this reference fe are $((:L2, )).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{CrouzeixRaviart}, conf::L2Conformity)
  get_face_dofs(reffe)
end

"""
    _get_dfaces_measure(p::Polytope{D}, d::Int)

Return the vector of `d`-volumes of the `d`-faces of `p`.
""" # TODO: Generalize
function _get_dfaces_measure(p::Polytope{D}, d::Int) where D
  @notimplementedif !is_simplex(p) || D>3 "only implemented for simplices of dim up to 3".
  measures = Float64[]
  dfaces_vertices = get_face_coordinates(p,d)
  for entity in dfaces_vertices
    n = length(entity)
    if n == 1  # The set containing one point has cardinal 1
       push!(measures, 1.0)
    elseif n == 2 # Length of an edge
      p1, p2 = entity
      push!(measures, norm(p2-p1))
    elseif n == 3 # Area of a simplex
      p1, p2, p3 = entity
      v1 = p2 - p1
      v2 = p3 - p1
      area = 0.5 * norm(cross(v1, v2))
      push!(measures, area)
    elseif n == 4 && d == 3 # Volume of a tetrahedron
      p1, p2, p3, p4 = entity
      v1 = p2 - p1
      v2 = p3 - p1
      v3 = p4 - p1
      volume = abs(dot(v1, cross(v2, v3))) / 6
      push!(measures, volume)
    end
  end
  return measures
end

