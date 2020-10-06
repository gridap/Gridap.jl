
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

function update_state_variables!(updater::Function,quad::CellQuadrature,f::CellField...)
  x = get_coordinates(quad)
  update_state_variables!(updater,x,f...)
end

function update_state_variables!(updater::Function,x::AbstractArray,f::CellField...)
  fx = map(i->evaluate(i,x),f)
  caches = map(array_cache,fx...)
  cache_x = array_cache(x)
  _update_state_variables!(updater,caches,fx,cache_x,x)
end

function update_state_variables!(quad::CellQuadrature,updater::Function,f::CellField...)
  msg =
  """
  The method
      update_state_variables!(quad::CellQuadrature,updater::Function,f::CellField...)
  has been removed. Use
      update_state_variables!(updater::Function,quad::CellQuadrature,f::CellField...)
  instead
  """
  error(msg)
end

@noinline function  _update_state_variables!(updater,caches,fx,cache_x,x)
  ncells = length(x)
  for cell in 1:ncells
    fxi = getitems!(caches,fx,cell)
    xi = getindex!(cache_x,x,cell)
    for q in 1:length(xi)
      fxiq = getitems(fxi,q)
      r = updater(fxiq...)
      need_to_update, states = Arrays._split(r...)
      if need_to_update
        _update_states!(fxi,q,states,Val{length(states)}())
      end
    end
  end
end

@inline function _update_states!(b,q,states,::Val{i}) where i
  _update_state!(b,q,states,Val{i}())
  _update_states!(b,q,states,Val{i-1}())
  nothing
end

@inline function _update_states!(b,q,states,::Val{0})
  nothing
end

@inline function _update_state!(b,q,states,::Val{i}) where i
  m = length(b)
  n = length(states)
  o = m-n
  b[i+o][q] = states[i]
  nothing
end
