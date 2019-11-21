

"""
    struct Table{T,D,P} <: AbstractVector{Vector{T}}
      data::D
      ptrs::P
    end

Type representing a list of lists (i.e., a table) in 
compressed format.
"""
struct Table{T,D,P} <: AbstractVector{Vector{T}}
  data::D
  ptrs::P

  @doc """
  """
  function Table(data::AbstractVector,ptrs::AbstractVector)
    T = eltype(data)
    D = typeof(data)
    P = typeof(ptrs)
    new{T,D,P}(data,ptrs)
  end
end

size(a::Table) = (length(a.ptrs)-1,)

IndexStyle(::Type{<:Table}) = IndexLinear()

function array_cache(a::Table)
  if length(a.ptrs) > 1
    pini = a.ptrs[1]
    pend = a.ptrs[2]
  else
    pend = a.ptrs[1]
    pini = pend
  end
  T = eltype(a.data)
  l = pend - pini
  r = zeros(T,l)
  CachedArray(r)
end

function getindex!(c,a::Table,i::Integer)
  pini = a.ptrs[i]
  l = a.ptrs[i+1] - pini
  setsize!(c,(l,))
  pini -= 1
  r = c.array
  for j in 1:l
     @inbounds r[j] = a.data[pini+j]
  end
  r
end

function getindex(a::Table,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

