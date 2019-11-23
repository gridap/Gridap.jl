

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

"""
"""
function Table(a::AbstractVector{<:AbstractVector})
  data, ptrs = generate_data_and_ptrs(a)
  Table(data,ptrs)
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

# Helper functions related with Tables

"""
Rewind the given vector.
"""
function rewind_ptrs!(ptrs::AbstractVector{<:Integer})
  @inbounds for i in (length(ptrs)-1):-1:1
    ptrs[i+1] = ptrs[i]
  end
  ptrs[1] = 1
end

"""
Given a vector of integers, mutate it from length state to pointer state.
"""
function length_to_ptrs!(ptrs::AbstractArray{<:Integer})
  ptrs[1] = 1
  @inbounds for i in 1:(length(ptrs)-1)
    ptrs[i+1] += ptrs[i]
  end
end

"""
Given a vector of vectors compute the corresponding data and and ptrs
"""
function generate_data_and_ptrs(vv::AbstractVector{<:AbstractVector{T}}) where T
  ptrs = Vector{T}(undef,length(vv)+1)
  _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  data = Vector{T}(undef,ndata)
  _generate_data_and_ptrs_fill_data!(data,vv)
  (data, ptrs)
end

function _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  for (i,v) in enumerate(vv)
    ptrs[i+1] = length(v)
  end
end

function _generate_data_and_ptrs_fill_data!(data,vv)
  k = 1
  for v in vv
    for vi in v
      data[k] = vi
      k += 1
    end
  end
end

"""
"""
function append_ptrs(pa::AbstractVector{T},pb::AbstractVector{T}) where T
  na = length(pa)-1
  nb = length(pb)-1
  p = zeros(T,(na+nb+1))
  _append_count!(p,pa,pb)
  length_to_ptrs!(p)
  p
end

function _append_count!(p,pa,pb)
  na = length(pa)-1
  for ca in 1:na
    p[ca+1] = pa[ca+1] - pa[ca]
  end
  nb = length(pb)-1
  for cb in 1:nb
    p[cb+1+na] = pb[cb+1] - pb[cb]
  end
end
