
"""
"""
struct DomainContribution <: GridapType
  dict::IdDict{Triangulation,AbstractArray}
end

DomainContribution() = DomainContribution(IdDict{Triangulation,AbstractArray}())

num_domains(a::DomainContribution) = length(a.dict)

get_domains(a::DomainContribution) = keys(a.dict)

function get_contribution(a::DomainContribution,trian::Triangulation)
  if haskey(a.dict,trian)
    return a.dict[trian]
  else
    for trian_a in get_domains(a) 
      if get_cell_node_ids(trian_a) == get_cell_node_ids(trian)
        if hasproperty(trian,:a)
          if hasproperty(trian.a,:subcells)
              if isequivtrian(trian,trian_a)
                return get_contribution(a,trian_a)
              end
          end
        elseif hasproperty(trian,:subfacets)
            if isequivtrian_subfacet(trian,trian_a)
              return get_contribution(a,trian_a)
            end
        end
      end
    end
  @unreachable """\n
  There is not contribution associated with the given mesh in this DomainContribution object.
  """
  end
end

function isequivtrian(trian1::AppendedTriangulation,trian2::AppendedTriangulation) 
  sc1 = trian1.a.subcells
  sc2 = trian2.a.subcells
  return (sc1.cell_to_points == sc2.cell_to_points && sc1.cell_to_bgcell == sc2.cell_to_bgcell && sc1.point_to_coords == sc2.point_to_coords && sc1.point_to_rcoords == sc2.point_to_rcoords )
end

function isequivtrian_subfacet(trian1,trian2) 
  sf1 = trian1.subfacets
  sf2 = trian2.subfacets
  return (sf1.facet_to_points == sf2.facet_to_points && sf1.facet_to_normal == sf2.facet_to_normal && sf1.point_to_coords == sf2.point_to_coords && sf1.point_to_rcoords == sf2.point_to_rcoords )
end

Base.getindex(a::DomainContribution,trian::Triangulation) = get_contribution(a,trian)

function add_contribution!(a::DomainContribution,trian::Triangulation,b::AbstractArray,op=+)

  S = eltype(b)
  if !(S<:AbstractMatrix || S<:AbstractVector || S<:Number || S<:ArrayBlock)
    @unreachable """\n
    You are trying to add a contribution with eltype $(S).
    Only cell-wise matrices, vectors, or numbers are accepted.

    Make sure that you are defining the terms in your weak form correclty.
    """
  end

  if length(a.dict) > 0
    T = eltype(first(values(a.dict)))
    if T <: AbstractMatrix || S<:(ArrayBlock{A,2} where A)
      @assert S<:AbstractMatrix || S<:(ArrayBlock{A,2} where A) """\n
      You are trying to add a contribution with eltype $(S) to a DomainContribution that
      stores cell-wise matrices.

      Make sure that you are defining the terms in your weak form correclty.
      """
    elseif T <: AbstractVector || S<:(ArrayBlock{A,1} where A)
      @assert S<:AbstractVector || S<:(ArrayBlock{A,1} where A) """\n
      You are trying to add a contribution with eltype $(S) to a DomainContribution that
      stores cell-wise vectors.

      Make sure that you are defining the terms in your weak form correclty.
      """
    elseif T <: Number
      @assert S<:Number """\n
      You are trying to add a contribution with eltype $(S) to a DomainContribution that
      stores cell-wise numbers.

      Make sure that you are defining the terms in your weak form correclty.
      """
    end
  end

  if haskey(a.dict,trian)
    a.dict[trian] = lazy_map(Broadcasting(op),a.dict[trian],b)
  else
    if op == +
     a.dict[trian] = b
    else
     a.dict[trian] = lazy_map(Broadcasting(op),b)
    end
  end
  a
end

Base.sum(a::DomainContribution)= sum(map(sum,values(a.dict)))

Base.copy(a::DomainContribution) = DomainContribution(copy(a.dict))

function (+)(a::DomainContribution,b::DomainContribution)
  c = copy(a)
  for (trian,array) in b.dict
    add_contribution!(c,trian,array)
  end
  c
end

function (-)(a::DomainContribution,b::DomainContribution)
  c = copy(a)
  for (trian,array) in b.dict
    add_contribution!(c,trian,array,-)
  end
  c
end

function (*)(a::Number,b::DomainContribution)
  c = DomainContribution()
  for (trian,array_old) in b.dict
    s = size(get_cell_map(trian))
    array_new = lazy_map(Broadcasting(*),Fill(a,s),array_old)
    add_contribution!(c,trian,array_new)
  end
  c
end

(*)(a::DomainContribution,b::Number) = b*a

function get_array(a::DomainContribution)
  @assert num_domains(a) == 1 """\n
  Method get_array(a::DomainContribution) can be called only
  when the DomainContribution object involves just one domain.
  """
  a.dict[first(keys(a.dict))]
end

abstract type Measure <: GridapType end

function integrate(f,b::Measure)
  @abstractmethod
end

function get_cell_quadrature(b::Measure)
  @abstractmethod
end

function get_cell_points(a::Measure)
  quad = get_cell_quadrature(a)
  return get_cell_points(quad)
end

function (*)(a::Integrand,b::Measure)
  integrate(a.object,b)
end

(*)(b::Measure,a::Integrand) = a*b

struct GenericMeasure <: Measure
  quad::CellQuadrature
end

Measure(q::CellQuadrature) = GenericMeasure(q)

Measure(args...;kwargs...) = GenericMeasure(CellQuadrature(args...;kwargs...))

get_cell_quadrature(a::GenericMeasure) = a.quad

function integrate(f,b::GenericMeasure)
  c = integrate(f,b.quad)
  cont = DomainContribution()
  add_contribution!(cont,b.quad.trian,c)
  cont
end

"""
  Composite Measure

  Measure such that the integration and target triangulations are different. 

  - ttrian: Target triangulation, where the domain contribution lives.
  - itrian: Integration triangulation, where the integration takes place.
  - quad  : CellQuadrature, defined in itrian
"""
struct CompositeMeasure{A<:Triangulation,B<:Triangulation,C<:CellQuadrature} <: Measure
  ttrian :: A
  itrian :: B
  quad   :: C
end

function Measure(ttrian::Triangulation,itrian::Triangulation,q::CellQuadrature)
  @check get_triangulation(q) === itrian
  return CompositeMeasure(ttrian,itrian,q)
end

Measure(ttrian::Triangulation,itrian::Triangulation,args...;kwargs...) = CompositeMeasure(ttrian,itrian,args...;kwargs...)

function CompositeMeasure(ttrian::Triangulation,itrian::Triangulation,args...;kwargs...)
  quad = CellQuadrature(itrian,args...;kwargs...)
  return CompositeMeasure(ttrian,itrian,quad)
end

get_cell_quadrature(a::CompositeMeasure) = a.quad

function integrate(f,b::CompositeMeasure)
  ic   = integrate(f,b.quad)
  cont = DomainContribution()
  tc   = move_contributions(ic,b.itrian,b.ttrian)
  add_contribution!(cont,b.ttrian,tc)
  return cont
end