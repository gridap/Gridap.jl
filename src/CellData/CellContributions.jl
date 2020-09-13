
struct LebesgueMeasure <: GridapType
  quad::CellQuadrature
  function LebesgueMeasure(quad::CellQuadrature)
    new(quad)
  end
end

function LebesgueMeasure(args...)
  quad = CellQuadrature(args...)
  LebesgueMeasure(quad)
end

function Base.:*(a::Integral,b::LebesgueMeasure)
  cell_values = a*b.quad
  cell_id = get_cell_id(get_coordinates(b.quad))
  CellContribution(cell_values,cell_id,b)
end

struct CellContribution{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  cell_values::A
  cell_id
  object
end

Base.size(a::CellContribution) = size(get_array(a))
Arrays.array_cache(a::CellContribution) = array_cache(get_array(a))
@inline Arrays.getindex!(cache,a::CellContribution,i::Integer...) = getindex!(cache,get_array(a),i...)
@inline Base.getindex(a::CellContribution,i::Integer...) = get_array(a)[i...]
Base.IndexStyle(::Type{CellContribution{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.sum(a::CellContribution) = sum(get_array(a))

Arrays.get_array(a::CellContribution) = a.cell_values
get_cell_id(a::CellContribution) = a.cell_id
get_object(a::CellContribution) = a.object

struct CollectionOfCellContribution <: GridapType
  dict::Dict{UInt,CellContribution}
end

function get_cell_contribution(a::CollectionOfCellContribution,object)
  id = objectid(object)
  a.dict[id]
end

function Base.:+(a::CellContribution,b::CellContribution)
  if get_object(a) === get_object(b)
    array = apply(bcast(+),get_array(a),get_array(b))
    CellContribution(array,get_cell_id(a),get_object(a))
  else
    dict = Dict{UInt,CellContribution}()
    dict[objectid(get_object(a))] = a
    dict[objectid(get_object(b))] = b
    CollectionOfCellContribution(dict)
  end
end

function Base.:+(a::CollectionOfCellContribution,b::CellContribution)
  id = objectid(get_object(b))
  if haskey(a.dict,id)
    c = a.dict[id]
    a.dict[id] = c + b
  else
    a.dict[id] = b
  end
  a
end

function Base.:+(a::CellContribution,b::CollectionOfCellContribution)
  b + a
end

