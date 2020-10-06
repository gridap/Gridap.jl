
"""
Array of vectors that combines an `AbstractArray{Array{T,M},N}` and another
array of masks with the same structure but `Bool` entries, and returns a
vector with the entries of the original array that are `true` in the mask array,
i.e., the filtered entries
"""
struct FilteredCellArray{T,M,N,L,V} <: AbstractArray{Array{T,M},N}
  cell_values::L
  cell_filters::V

  @doc """
  """
  function FilteredCellArray(
    cell_values::AbstractArray{<:AbstractArray},
    cell_filters::AbstractArray{<:AbstractArray})

    @assert size(cell_values) == size(cell_filters) "Global arrays mismatch"

    T = eltype(eltype(cell_values))

    L = typeof(cell_values)
    V = typeof(cell_filters)

    # M = ndims(eltype(L))
    M = 1 # The result is always a vector with the entries
    N = ndims(L)

    O = ndims(V)

    @assert N == O "Local arrays dim mismatch"

    @assert (eltype(eltype(cell_filters))<: Bool) "Filters are not Booleans"

    new{T,M,N,L,V}(cell_values,cell_filters)

  end
end

size(a::FilteredCellArray) = size(a.cell_values)

function IndexStyle(::Type{<:FilteredCellArray{T,M,N,L}}) where {T,M,N,L}
  IndexStyle(L)
end

function array_cache(a::FilteredCellArray)
  vals = testitem(a.cell_values)
  T = eltype(eltype(a.cell_values))
  r = zeros(T,length(vals))
  c = CachedArray(r)
  cv = array_cache(a.cell_values)
  cb = array_cache(a.cell_filters)
  (cv,cb,c)
end

function getindex!(cache,a::FilteredCellArray,i::Integer...)
  (cv,cb,c) = cache
  vals = getindex!(cv,a.cell_values,i...)
  filters = getindex!(cb,a.cell_filters,i...)
  @assert size(vals) == size(filters) "Local arrays mismatch"
  setsize!(c,(sum(filters),))
  r = c.array
  i = 0
  for (val,filter) in zip(vals,filters)
    if filter
      i += 1
      r[i] = val
    end
  end
  r
end

function getindex(a::FilteredCellArray,i::Integer...)
  cache = array_cache(a)
  getindex!(cache,a,i...)
end

struct FilterKernel <: Kernel end

function kernel_return_type(f::FilterKernel,x...)
  typeof(kernel_testitem(f,x...))
end

function kernel_cache(k::FilterKernel,f,a)
  # vals = testitem(a)
  vals = a
  T = eltype(eltype(a))
  r = zeros(T,length(vals))
  c = CachedArray(r)
end

function lazy_map_kernel!(cache,k::FilterKernel,f,a)
  c = cache
  vals = a
  filters = f
  @assert size(vals) == size(filters) "Local arrays mismatch"
  setsize!(c,(sum(filters),))
  r = c.array
  i = 0
  for (val,filter) in zip(vals,filters)
    if filter
      i += 1
      r[i] = val
    end
  end
  r
end
