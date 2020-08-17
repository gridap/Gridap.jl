
"""
"""
function QPointCellField(value::Number,cell_map::AbstractArray{<:Field},quad::CellQuadrature)
  q = get_coordinates(quad)
  v_q = [ fill(value,size(qi)) for qi in q ]
  array = ArrayOfEvaluatedFields(v_q,q)
  GenericCellField(array, cell_map)
end

"""
"""
function CellField(value::Number,cell_map::AbstractArray{<:Field},quad::CellQuadrature)
  QPointCellField(value,cell_map,quad)
end

struct EvaluatedField{A<:AbstractArray,P<:AbstractArray} <:Field
  array::A
  points::P
end

function field_cache(f::EvaluatedField,x)
  nothing
end

function evaluate_field!(cache,f::EvaluatedField,x)
  @assert length(x) == length(f.array)
  @assert x === f.points || x == f.points
  f.array
end

struct ArrayOfEvaluatedFields{T,N,A,B} <: AbstractArray{EvaluatedField{T},N}
  array::A
  points::B
  function ArrayOfEvaluatedFields(array::AbstractArray{T,N},points::AbstractArray) where {T<:AbstractArray,N}
    A = typeof(array)
    B = typeof(points)
    new{T,N,A,B}(array,points)
  end
end

Base.size(a::ArrayOfEvaluatedFields) = size(a.array)

Base.IndexStyle(::Type{<:ArrayOfEvaluatedFields{T,N,A}}) where {T,N,A} = IndexStyle(A)

@inline function Base.getindex(a::ArrayOfEvaluatedFields,i::Integer)
  EvaluatedField(a.array[i],a.points[i])
end

@inline function Base.getindex(a::ArrayOfEvaluatedFields{T,N},i::Vararg{Int,N}) where {T,N}
  EvaluatedField(a.array[i...],a.points[i...])
end

function array_cache(a::ArrayOfEvaluatedFields)
  ca = array_cache(a.array)
  cp = array_cache(a.points)
  (ca,cp)
end

@inline function getindex!(cache,a::ArrayOfEvaluatedFields,i...) 
  ca, cp = cache
  array = getindex!(ca,a.array,i...)
  points = getindex!(ca,a.points,i...)
  EvaluatedField(array,points)
end

function evaluate_field_array(a::ArrayOfEvaluatedFields,x::AbstractArray)
  @assert a.points === x || a.points == x
  @assert length(a) == length(x)
  a.array
end


