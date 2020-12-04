
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
    @unreachable """\n
    There is not contribution associated with the given mesh in this DomainContribution object.
    """
  end
end

Base.getindex(a::DomainContribution,trian::Triangulation) = get_contribution(a,trian)

function add_contribution!(a::DomainContribution,trian::Triangulation,b::AbstractArray,op=+)

  S = eltype(b)
  if !(S<:AbstractMatrix || S<:AbstractVector || S<:Number)
    @unreachable """\n
    You are trying to add a contribution with eltype $(S).
    Only cell-wise matrices, vectors, or numbers are accepted.

    Make sure that you are defining the terms in your weak form correclty.
    """
  end

  if length(a.dict) > 0
    T = eltype(first(values(a.dict)))
    if T <: AbstractMatrix
      @assert S<:AbstractMatrix """\n
      You are trying to add a contribution with eltype $(S) to a DomainContribution that
      stores cell-wise matrices.

      Make sure that you are defining the terms in your weak form correclty.
      """
    elseif T <: AbstractVector
      @assert S<:AbstractVector """\n
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

function get_array(a::DomainContribution)
  @assert num_domains(a) == 1 """\n
  Method get_array(a::DomainContribution) can be called only
  when the DomainContribution object involves just one domain.
  """
  a.dict[first(keys(a.dict))]
end

struct LebesgueMeasure <: GridapType
  quad::CellQuadrature
end

get_cell_points(a::LebesgueMeasure) = get_cell_points(a.quad)

LebesgueMeasure(args...) = LebesgueMeasure(CellQuadrature(args...))

function integrate(f,b::LebesgueMeasure)
  c = integrate(f,b.quad)
  cont = DomainContribution()
  add_contribution!(cont,b.quad.trian,c)
  cont
end

function (*)(a::Integrand,b::LebesgueMeasure)
  integrate(a.object,b)
end

(*)(b::LebesgueMeasure,a::Integrand) = a*b


