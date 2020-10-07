
"""
    struct CompressedArray{T,N,A,P} <: AbstractArray{T,N}
      values::A
      ptrs::P
    end
Type representing an array with a reduced set of values.
The array is represented by a short array of values, namely
the field `values`, and a large array of indices, namely the
field `ptrs`. The `i`-th component of the resulting array is
defined as `values[ptrs[i]]`. The type parameters `A`, and `P`
are restricted to be array types by the inner constructor of this `struct`.
"""
struct CompressedArray{T,N,A,P} <: AbstractArray{T,N}
  values::A
  ptrs::P
  @doc """
      CompressedArray(values::AbstractArray,ptrs::AbstractArray)

  Creates a `CompressedArray` object by the given arrays of `values` and
  `ptrs`.
  """
  function CompressedArray(values::AbstractArray,ptrs::AbstractArray)
    A = typeof(values)
    P = typeof(ptrs)
    T = eltype(values)
    N = ndims(ptrs)
    new{T,N,A,P}(values,ptrs)
  end
end

function testitem(a::CompressedArray)
  if length(a.ptrs) == 0
    testitem(a.values)
  else
    a.values[first(a.ptrs)]
  end
end

size(a::CompressedArray) = size(a.ptrs)

@propagate_inbounds function getindex(a::CompressedArray,i::Integer)
  j = a.ptrs[i]
  a.values[j]
end

@propagate_inbounds function getindex(a::CompressedArray,i::Integer...)
  j = a.ptrs[i...]
  a.values[j]
end

function IndexStyle(a::Type{CompressedArray{T,N,A,P}}) where {T,N,A,P}
  IndexStyle(P)
end

function lazy_map(::typeof(evaluate),f::Fill,g1::CompressedArray,g::CompressedArray...)
  if all( ( gi.ptrs === g1.ptrs for gi in g ) ) || all( ( gi.ptrs == g1.ptrs for gi in g ) )
    _lazy_map_fill_compressed(f,g1,g...)
  else
    return LazyArray(f,g1,g...)
  end
end

function lazy_map(::typeof(evaluate),g1::CompressedArray,g::CompressedArray...)
  if all( ( gi.ptrs === g1.ptrs for gi in g ) ) || all( ( gi.ptrs == g1.ptrs for gi in g ) )
    _lazy_map_compressed(g1,g...)
  else
    return LazyArray(g1,g...)
  end
end

function lazy_map(::typeof(evaluate),g1::CompressedArray,g::Fill...)
  f = _fill_to_compressed(g1,g)
  _lazy_map_compressed(g1,f...)
end

function lazy_map(::typeof(evaluate),f::Fill,g1::CompressedArray,g::Fill...)
  h = _fill_to_compressed(g1,g)
  _lazy_map_fill_compressed(f,g1,h...)
end

function lazy_map(::typeof(evaluate),f::Fill,g1::CompressedArray)
  _lazy_map_fill_compressed(f,g1)
end

function lazy_map(::typeof(evaluate),::Type{T},f::Fill,g1::CompressedArray,g::CompressedArray...) where T
  if all( ( gi.ptrs === g1.ptrs for gi in g ) ) || all( ( gi.ptrs == g1.ptrs for gi in g ) )
    _lazy_map_fill_compressed(f,g1,g...)
  else
    return LazyArray(T,f,g1,g...)
  end
end

function lazy_map(::typeof(evaluate),::Type{T},g1::CompressedArray,g::CompressedArray...) where T
  if all( ( gi.ptrs === g1.ptrs for gi in g ) ) || all( ( gi.ptrs == g1.ptrs for gi in g ) )
    _lazy_map_compressed(g1,g...)
  else
    return LazyArray(T,g1,g...)
  end
end

function lazy_map(::typeof(evaluate),::Type{T},g1::CompressedArray,g::Fill...) where T
  f = _fill_to_compressed(g1,g)
  _lazy_map_compressed(g1,f...)
end

function lazy_map(::typeof(evaluate),::Type{T},f::Fill,g1::CompressedArray,g::Fill...) where T
  h = _fill_to_compressed(g1,g)
  _lazy_map_fill_compressed(f,g1,h...)
end

function lazy_map(::typeof(evaluate),::Type{T},f::Fill,g1::CompressedArray) where T
  _lazy_map_fill_compressed(f,g1)
end

function _fill_to_compressed(g1,g)
  ptrs = g1.ptrs
  l = length(g1.values)
  f = ( CompressedArray(fill(gi.value,l),ptrs) for gi in g )
  f
end

function _lazy_map_fill_compressed(f,g1,g...)
  k = f.value
  ptrs = g1.ptrs
  vals = _getvalues(g1,g...)
  vk = lazy_map(k,vals...)
  CompressedArray(collect(vk),ptrs)
end

function _lazy_map_compressed(g1,g...)
  ptrs = g1.ptrs
  vals = _getvalues(g...)
  vk = lazy_map(evaluate,g1.values,vals...)
  CompressedArray(collect(vk),ptrs)
end

function _getvalues(a,b...)
  va = a.values
  vb = _getvalues(b...)
  (va,vb...)
end

function _getvalues(a)
  va = a.values
  (va,)
end
