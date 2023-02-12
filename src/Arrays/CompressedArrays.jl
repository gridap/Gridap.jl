
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

@propagate_inbounds function getindex(a::CompressedArray{T,N},i::Vararg{Integer,N}) where {T,N}
  j = a.ptrs[i...]
  a.values[j]
end

function IndexStyle(a::Type{CompressedArray{T,N,A,P}}) where {T,N,A,P}
  IndexStyle(P)
end

function lazy_map(::typeof(evaluate),::Type{T},g::CompressedArray...) where T
  if _have_same_ptrs(g)
    _lazy_map_compressed(g...)
  else
    LazyArray(T,g...)
  end
end

function lazy_map(::typeof(evaluate),::Type{T},g::Union{CompressedArray,Fill}...) where T
  g_compressed = _find_compressed_ones(g)
  if _have_same_ptrs(g_compressed)
    g1 = first(g_compressed)
    g_all_compressed = map(gi->_compress(gi,g1),g)
    _lazy_map_compressed(g_all_compressed...)
  else
    LazyArray(T,g...)
  end
end

function _find_compressed_ones(g)
  g_compressed = ( gi for gi in g if isa(gi,CompressedArray) )
  g_compressed
end

function _lazy_map_compressed(g::CompressedArray...)
  vals = map(evaluate, map(gi->gi.values,g)...)
  ptrs = first(g).ptrs
  CompressedArray(vals,ptrs)
end

function _have_same_ptrs(g)
  g1 = first(g)
  all(map( gi -> gi.ptrs === g1.ptrs ||  gi.ptrs == g1.ptrs, g))
end

function _compress(a::CompressedArray,b::CompressedArray)
  @check _have_same_ptrs((a,b))
  a
end

function _compress(a::Fill,b::CompressedArray)
  vals = fill(a.value,length(b.values))
  CompressedArray(vals,b.ptrs)
end

function same_branch(a::CompressedArray,b::CompressedArray)
  a.ptrs === b.ptrs && a.values == b.values
end

function Base.unique(a::CompressedArray)
  unique(a.values)
end
