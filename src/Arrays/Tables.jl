

"""
     struct Table{T,Vd<:AbstractVector{T},Vp<:AbstractVector} <: AbstractVector{Vector{T}}
        data::Vd
        ptrs::Vp
     end

Type representing a list of lists (i.e., a table) in
compressed format.
"""
struct Table{T,Vd<:AbstractVector{T},Vp<:AbstractVector} <: AbstractVector{Vector{T}}
  data::Vd
  ptrs::Vp
  function Table(data::AbstractVector,ptrs::AbstractVector)
    new{eltype(data),typeof(data),typeof(ptrs)}(data,ptrs)
  end
end

"""
    Table(a::AbstractArray{<:AbstractArray})

Build a table from a vector of vectors. If the inputs are
multidimensional arrays instead of vectors, they are flattened.
"""
function Table(a::AbstractArray{<:AbstractArray})
  data, ptrs = generate_data_and_ptrs(a)
  Table(data,ptrs)
end

function Table(a::Table)
  a
end

function Base.convert(::Type{Table{T,Vd,Vp}},table::Table{Ta,Vda,Vpa}) where {T,Vd,Vp,Ta,Vda,Vpa}
  data = convert(Vd,table.data)
  ptrs = convert(Vp,table.ptrs)
  Table(data,ptrs)
end

function Base.convert(::Type{Table{T,Vd,Vp}},table::Table{T,Vd,Vp}) where {T,Vd,Vp}
  table
end

function Base.view(a::Table,i::Integer)
  pini = a.ptrs[i]
  pend = a.ptrs[i+1]-1
  return view(a.data,pini:pend)
end

function Base.view(a::Table,ids::UnitRange{<:Integer})
  data_range = a.ptrs[ids.start]:a.ptrs[ids.stop+1]-1
  ptrs_range = ids.start:ids.stop+1
  offset = a.ptrs[ids.start]-1
  ptrs = lazy_map(p -> p - offset, view(a.ptrs,ptrs_range))
  return Table(view(a.data,data_range),ptrs)
end

"""
"""
function identity_table(::Type{T},::Type{P},l::Integer) where {T,P}
  data = Vector{T}(1:l)
  ptrs = Vector{P}(1:l+1)
  Table(data,ptrs)
end

"""
    empty_table(::Type{T},::Type{P}, l::Integer) where {T,P}
    empty_table(l::Integer)
"""
function empty_table(::Type{T},::Type{P}, l::Integer) where {T,P}
  data = T[]
  ptrs = ones(P,l+1)
  Table(data,ptrs)
end

empty_table(l::Integer) = empty_table(Int,Int32,l)


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
  r = Vector{T}(undef,l)
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

function Base.getindex(a::Table,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function Base.getindex(a::Table,i::UnitRange)
  r = a.ptrs[i.start]:(a.ptrs[i.stop+1]-1)
  data = a.data[r]
  r = i.start:(i.stop+1)
  ptrs = a.ptrs[r]
  o = ptrs[1]-1
  ptrs .-= o
  Table(data,ptrs)
end

function Base.getindex(a::Table,ids::AbstractVector{<:Integer})
  ptrs = similar(a.ptrs,eltype(a.ptrs),length(ids)+1)
  for (i,id) in enumerate(ids)
    ptrs[i+1] = a.ptrs[id+1]-a.ptrs[id]
  end
  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  data = similar(a.data,eltype(a.data),ndata)
  for (i,id) in enumerate(ids)
    n = a.ptrs[id+1]-a.ptrs[id]
    p1 = ptrs[i]-1
    p2 = a.ptrs[id]-1
    for j in 1:n
      data[p1+j] = a.data[p2+j]
    end
  end
  Table(data,ptrs)
end

# Helper functions related with Tables

"""
    data, ptrs = generate_data_and_ptrs(vv)

Given a vector of vectors, compress it and return the corresponding data and and ptrs
"""
function generate_data_and_ptrs(vv::AbstractArray{<:AbstractArray{T}}) where T
  ptrs = Vector{Int32}(undef,length(vv)+1)
  _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  data = Vector{T}(undef,ndata)
  _generate_data_and_ptrs_fill_data!(data,vv)
  (data, ptrs)
end

function _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  c = array_cache(vv)
  k = 1
  for i in eachindex(vv)
    v = getindex!(c,vv,i)
    ptrs[k+1] = length(v)
    k += 1
  end
end

function _generate_data_and_ptrs_fill_data!(data,vv)
  c = array_cache(vv)
  k = 1
  for i in eachindex(vv)
    v = getindex!(c,vv,i)
    for vi in v
      data[k] = vi
      k += 1
    end
  end
end

"""
    append_ptrs(pa,pb)

Concatenate two vectors of pointers in a new vector.
"""
function append_ptrs(pa::AbstractVector{T},pb::AbstractVector{T}) where T
  p = copy(pa)
  append_ptrs!(p,pb)
end

"""
    append_ptrs!(pa,pb)

Similar to [`append_ptrs`](@ref), but appends `pb` at the end of `pa`, in place
in `pa`.
"""
function append_ptrs!(pa::AbstractVector{T},pb::AbstractVector{T}) where T
  na = length(pa)-1
  nb = length(pb)-1
  _append_grow!(pa,nb,zero(T))
  _append_count!(pa,pb,na,nb)
  rewind_ptrs!(pa)
  length_to_ptrs!(pa)
  pa
end

function _append_grow!(pa,nb,z)
  for i in 1:nb
    push!(pa,z)
  end
end

function _append_count!(pa,pb,na,nb)
  for ca in 1:na
    pa[ca] = pa[ca+1] - pa[ca]
  end
  for cb in 1:nb
    pa[cb+na] = pb[cb+1] - pb[cb]
  end
end

"""
"""
const UNSET = 0

"""
"""
function find_inverse_index_map(a_to_b, nb=maximum(a_to_b))
  T = eltype(a_to_b)
  b_to_a = fill(T(UNSET),nb)
  find_inverse_index_map!(b_to_a, a_to_b)
  b_to_a
end

"""
"""
function find_inverse_index_map!(b_to_a, a_to_b)
  for (a,b) in enumerate(a_to_b)
    if b != UNSET
      @inbounds b_to_a[b] = a
    end
  end
end

"""
"""
function append_tables_globally(tables::Table{T,Vd,Vp}...) where {T,Vd,Vp}
  first_table, = tables
  data = copy(first_table.data)
  ptrs = copy(first_table.ptrs)
  for (i,table) in enumerate(tables)
    if  i != 1
      append!(data,table.data)
      append_ptrs!(ptrs,table.ptrs)
    end
  end
  Table(data,ptrs)
end

function append_tables_globally()
  @unreachable "At least one table has to be provided"
end

"""
"""
get_ptrs_eltype(::Table{T,Vd,Vp}) where {T,Vd,Vp} = eltype(Vp)
get_ptrs_eltype(::Type{Table{T,Vd,Vp}}) where {T,Vd,Vp} = eltype(Vp)

"""
"""
get_data_eltype(::Table{T,Vd,Vp}) where {T,Vd,Vp} = T
get_data_eltype(::Type{Table{T,Vd,Vp}}) where {T,Vd,Vp} = T

"""
    append_tables_locally(tables::Table...)
"""
function append_tables_locally(tables::Table...)
  n = length(tables)
  offsets = tfill(0,Val{n}())
  append_tables_locally(offsets,tables)
end

"""
"""
function append_tables_locally(offsets::NTuple, tables::NTuple)

  first_table, = tables

  @check all( map(length,tables) .== length(first_table) ) "All tables must have the same length"
  ndata = sum( (length(table.data) for table in tables) )

  T = get_data_eltype(first_table)
  P = get_ptrs_eltype(first_table)

  ptrs = zeros(P,length(first_table.ptrs))
  data = zeros(T,ndata)

  for table in tables
    _append_tables_locally_count!(ptrs,table)
  end

  length_to_ptrs!(ptrs)

  for (offset,table) in zip(offsets,tables)
    _append_tables_locally_fill!(data,ptrs,offset,table)
  end

  rewind_ptrs!(ptrs)

  Table(data,ptrs)

end

function append_tables_locally(offsets::Tuple{}, tables::Tuple{})
  @unreachable "At least one table has to be provided"
end

function  _append_tables_locally_count!(ptrs,table)
  table_ptrs = table.ptrs
  n = length(table_ptrs)-1
  for i in 1:n
    ptrs[i+1] += table_ptrs[i+1]- table_ptrs[i]
  end
end

function _append_tables_locally_fill!(data,ptrs,offset,table)
  table_data = table.data
  table_ptrs = table.ptrs
  n = length(table_ptrs)-1
  for i in 1:n
    for j in table_ptrs[i]:(table_ptrs[i+1]-1)
      data[ptrs[i]] = table_data[j] + offset
      ptrs[i] += 1
    end
  end
end

"""
    collect1d(a)

Equivalent to

    [a[i] for in 1:length(a)]
"""
collect1d(a) = [a[i] for i in 1:length(a)]

function lazy_map(::typeof(getindex),a::Table,b::AbstractArray{<:Integer})
  LocalItemFromTable(a,b)
end

"""
    get_local_item(a::Table,li::Integer)

View in the `li`ᵗʰ column of `a` (the `li`ᵗʰ items in each list/row of `a`).
"""
function get_local_item(a::Table,li::Integer)
  LocalItemFromTable(a,Fill(li,length(a)))
end

struct LocalItemFromTable{T,Vd,Vp,A} <: AbstractVector{T}
  a_to_lb_to_b::Table{T,Vd,Vp}
  a_to_lb::A
end

Base.size(m::LocalItemFromTable) = size(m.a_to_lb_to_b)

Base.IndexStyle(::Type{<:LocalItemFromTable}) = IndexLinear()

@propagate_inbounds function Base.getindex(m::LocalItemFromTable, a::Integer)
  p = m.a_to_lb_to_b.ptrs[a]-1
  m.a_to_lb_to_b.data[p+m.a_to_lb[a]]
end

"""
    find_local_index(a_to_b, b_to_la_to_a)
"""
function find_local_index(a_to_b, b_to_la_to_a)
  @notimplemented "find_local_index only implemented for table"
end

function find_local_index(a_to_b, b_to_la_to_a::Table)
  a_to_la = LocalIndexFromTable(a_to_b, b_to_la_to_a)
  a_to_la
end

struct LocalIndexFromTable{T,Vd,Vp,V<:AbstractVector} <: AbstractVector{T}
  a_to_b::V
  b_to_la_to_a::Table{T,Vd,Vp}
end

Base.size(m::LocalIndexFromTable) = size(m.a_to_b)

Base.IndexStyle(::Type{<:LocalIndexFromTable}) = IndexStyle(Table)

@propagate_inbounds function Base.getindex(m::LocalIndexFromTable{T}, a::Integer) where T
  b = m.a_to_b[a]
  pini = m.b_to_la_to_a.ptrs[b]
  pend = m.b_to_la_to_a.ptrs[b+1]-1
  la = zero(T)
  for (la,p) in enumerate(pini:pend)
    if a == m.b_to_la_to_a.data[p]
      return T(la)
    end
  end
  return T(UNSET)
end

"""
    flatten_partition(a_to_bs::Table,nb::Integer)
    flatten_partition(a_to_bs::Table)
"""
function flatten_partition(a_to_bs::Table,nb::Integer=maximum(a_to_bs.data))
  T = eltype(eltype(a_to_bs))
  b_to_a = zeros(T,nb)
  flatten_partition!(b_to_a,a_to_bs)
  b_to_a
end

function  flatten_partition!(b_to_a,a_to_bs::Table)
  for a in 1:length(a_to_bs)
    pini = a_to_bs.ptrs[a]
    pend = a_to_bs.ptrs[a+1]-1
    for p in pini:pend
      b = a_to_bs.data[p]
      b_to_a[b] = a
    end
  end
end

function to_dict(table::Table)
  dict = Dict{Symbol,Any}()
  dict[:data] = table.data
  dict[:ptrs] = table.ptrs
  dict
end

function from_dict(::Type{Table{T,Vd,Vp}}, dict::Dict{Symbol,Any}) where {T,Vd,Vp}
  data::Vd = dict[:data]
  ptrs::Vp = dict[:ptrs]
  Table(data,ptrs)
end

function Base.copy(a::Table)
  Table(copy(a.data),copy(a.ptrs))
end
