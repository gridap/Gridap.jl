

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
"""
get_ptrs_eltype(::Table{T,Vd,Vp}) where {T,Vd,Vp} = eltype(Vp)
get_ptrs_eltype(::Type{Table{T,Vd,Vp}}) where {T,Vd,Vp} = eltype(Vp)

"""
"""
get_data_eltype(::Table{T,Vd,Vp}) where {T,Vd,Vp} = T
get_data_eltype(::Type{Table{T,Vd,Vp}}) where {T,Vd,Vp} = T

"""
    Table(a::AbstractArray{<:AbstractArray})

Build a table from a vector of vectors. If the inputs are
multidimensional arrays instead of vectors, they are flattened.
"""
function Table(a::AbstractArray{<:AbstractArray})
  data, ptrs = generate_data_and_ptrs(a)
  Table(data,ptrs)
end

Table(a::Table) = a

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

"""
    datarange(a::Table,i::Integer)
    datarange(a::Table,ids::UnitRange{<:Integer})

Given a `Table` and an index or a range of indices, return the
corresponding range of indices in the underlying data array.
Similar to `nzrange` for sparse matrices, it allows for convenient 
iteration over a table:

```julia
t = Table([[4,7],[8],[9,20,1]])
for i in eachindex(t)
  for k in datarange(t,i)
    val = t.data[k]
    # stuff ... 
  end
end
```

"""
@inline function datarange(a::Table,i::Integer)
  pini = a.ptrs[i]
  pend = a.ptrs[i+1]-1
  return pini:pend
end

@inline function datarange(a::Table,ids::UnitRange{<:Integer})
  pini = a.ptrs[ids.start]
  pend = a.ptrs[ids.stop+1]-1
  return pini:pend
end

"""
    dataview(a::Table, i::Integer)
    dataview(a::Table, ids::UnitRange{<:Integer})

Given a `Table` and an index or a range of indices, return a view
of the corresponding entries in the underlying data array.

Equivalent to `view(a.data, datarange(a,i))`.
"""
dataview(a::Table, ids) = view(a.data, datarange(a,ids))

"""
    dataiterator(a::Table)

Iterate over the entries of `a` returning the triplets `(i,j,v)` where 

- `i` is the outer index, 
- `j` is the local inner index, and
- `v` is the value `a[i][j]`.

## Example

```jldoctest
julia> t = Table([[4.,7.],[8.],[9.,2.,1.]])

julia> x = collect(dataiterator(t))

julia> x
6-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1, 4.)
 (1, 2, 7.)
 (2, 1, 8.)
 (3, 1, 9.)
 (3, 2, 2.)
 (3, 3, 1.)
```
"""
dataiterator(a::Table) = TableDataIterator(a)

struct TableDataIterator{T,Vd,Vp}
  t::Table{T,Vd,Vp}
end

Base.length(a::TableDataIterator) = length(a.t.data)
Base.size(a::TableDataIterator) = (length(a.t.data),)
Base.eltype(::Type{<:TableDataIterator{T}}) where T = Tuple{Int,Int,T}
Base.IteratorSize(::Type{<:TableDataIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:TableDataIterator}) = Base.HasEltype()

function Base.iterate(a::TableDataIterator)
  isempty(a.t) && (return nothing)
  i, j, v = 1, 1 , a.t.data[1]
  return (i,j,v), (i,j)
end

function Base.iterate(a::TableDataIterator, state)
  i, j = state
  data = a.t.data
  ptrs = a.t.ptrs
  n = length(ptrs)
  while (i < n) && (j+1 > ptrs[i+1] - ptrs[i])
    i += 1
    j = 0
  end
  (i == n) && (return nothing)
  v = data[ptrs[i] + j]
  j += 1
  return (i,j,v), (i,j)
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

function generate_ptrs(vv::AbstractArray{<:AbstractArray{T}}) where T
  ptrs = Vector{Int32}(undef,length(vv)+1)
  _generate_data_and_ptrs_fill_ptrs!(ptrs,vv)
  length_to_ptrs!(ptrs)
  ptrs
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

Append two vectors of pointers.
"""
function append_ptrs(pa::AbstractVector{T},pb::AbstractVector{T}) where T
  p = copy(pa)
  append_ptrs!(p,pb)
end

"""
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
    find_inverse_index_map(a_to_b[, nb=maximum(a_to_b)])
    find_inverse_index_map!(b_to_a, a_to_b)

Given a vector of indices `a_to_b`, returns the inverse index map `b_to_a`.
"""
function find_inverse_index_map(a_to_b, nb=maximum(a_to_b))
  T = eltype(a_to_b)
  b_to_a = fill(T(UNSET),nb)
  find_inverse_index_map!(b_to_a, a_to_b)
  b_to_a
end

function find_inverse_index_map!(b_to_a, a_to_b)
  for (a,b) in enumerate(a_to_b)
    if b != UNSET
      @inbounds b_to_a[b] = a
    end
  end
end

"""
    inverse_table(a_to_lb_to_b::Table [, nb=maximum(a_to_lb_to_b.data)])
    inverse_table(a_to_b::AbstractVector [, nb=maximum(a_to_b)])

Returns the inverse of the input `Table` or non-injective array of integers, as a `Table`.
"""
function inverse_table(
  a_to_lb_to_b::Table{T}, nb = maximum(a_to_lb_to_b.data,init=zero(T))
) where T
  data, ptrs = inverse_table(
    a_to_lb_to_b.data, a_to_lb_to_b.ptrs, nb
  )
  a_to_lb_to_b = Table(data,ptrs)
  return a_to_lb_to_b
end

function inverse_table(
  a_to_b::AbstractVector{T}, nb = maximum(a_to_b,init=zero(T))
) where T
  na = length(a_to_b)
  data, ptrs = inverse_table(
    a_to_b, 1:(na+1), nb
  )
  a_to_lb_to_b = Table(data,ptrs)
  return a_to_lb_to_b
end

function inverse_table(
  a_to_lb_to_b_data::AbstractVector{L},
  a_to_lb_to_b_ptrs::AbstractVector{P},
  nb::Integer
) where {L,P}

  o = one(P)
  ptrs = zeros(P,nb+1)
  @inbounds for b in a_to_lb_to_b_data
    ptrs[b+1] += o
  end
  length_to_ptrs!(ptrs)

  na = length(a_to_lb_to_b_ptrs)-1
  data = Vector{L}(undef,ptrs[end]-1)
  @inbounds for a in 1:na
    s = a_to_lb_to_b_ptrs[a]
    e = a_to_lb_to_b_ptrs[a+1] - o
    @inbounds for p in s:e
      b = a_to_lb_to_b_data[p]
      if b != UNSET
        data[ptrs[b]] = a
        ptrs[b] += o
      end
    end
  end
  rewind_ptrs!(ptrs)

  return data, ptrs
end

"""
"""
function append_tables_globally(
  first_table::Table{T,Vd,Vp}, tables::Table{T,Vd,Vp}...
) where {T,Vd,Vp}
  data = copy(first_table.data)
  ptrs = copy(first_table.ptrs)
  for table in tables
    append!(data,table.data)
    append_ptrs!(ptrs,table.ptrs)
  end
  Table(data,ptrs)
end

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
  @check length(offsets) == length(tables) !== 0 "Offsets and tables must have the same length"
  first_table, = tables
  @check all(t -> length(t) == length(first_table), tables) "All tables must have the same length"
  ndata = sum(t -> length(t.data), tables)

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
    remove_empty_entries!(table::Table)

Given a `Table`, remove the entries that are empty by modifying its `ptrs` in-place.
"""
function remove_empty_entries!(table::Table)
  ptrs = table.ptrs
  
  i = 1
  n = length(table)
  while i <= n
    if ptrs[i] == ptrs[i+1]
      deleteat!(ptrs,i+1)
      n -= 1
    else
      i += 1
    end
  end

  return table
end

"""
    collect1d(a)

Equivalent to

    [a[i] for in 1:length(a)]
"""
collect1d(a) = [a[i] for i in 1:length(a)]
collect1d(a::Vector) = a

function lazy_map(::typeof(getindex),a::Table,b::AbstractArray{<:Integer})
  LocalItemFromTable(a,b)
end

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
    find_local_index(a_to_b, b_to_la_to_a) -> a_to_la
    find_local_index(c_to_a, c_to_b, b_to_la_to_a) -> c_to_lc
"""
@inline function find_local_index(args...) 
  find_local_index(GridapLocalInt, args...)
end

function find_local_index(::Type{T}, a_to_b, b_to_la_to_a::Table) where T
  a_to_la = LocalIndexFromTable(T, a_to_b, b_to_la_to_a)
  a_to_la
end

struct LocalIndexFromTable{T,TT,Vd,Vp,V<:AbstractVector} <: AbstractVector{T}
  a_to_b::V
  b_to_la_to_a::Table{TT,Vd,Vp}

  function LocalIndexFromTable(
    ::Type{T}, a_to_b::AbstractVector, b_to_la_to_a::Table{TT,Vd,Vp}
  ) where {T,TT,Vd,Vp}
    V = typeof(a_to_b)
    new{T,TT,Vd,Vp,V}(a_to_b,b_to_la_to_a)
  end
end

LocalIndexFromTable(a_to_b, b_to_la_to_a) = LocalIndexFromTable(GridapLocalInt, a_to_b, b_to_la_to_a)

Base.size(m::LocalIndexFromTable) = size(m.a_to_b)

Base.IndexStyle(::Type{<:LocalIndexFromTable}) = IndexStyle(Table)

@propagate_inbounds function Base.getindex(m::LocalIndexFromTable{T}, a::Integer) where T
  b = m.a_to_b[a]
  for (la,p) in enumerate(datarange(m.b_to_la_to_a,b))
    if a == m.b_to_la_to_a.data[p]
      return T(la)
    end
  end
  return T(UNSET)
end

function find_local_index(
  ::Type{T}, c_to_a :: AbstractVector, c_to_b :: AbstractVector, b_to_la_to_a :: Table
) where T
  c_to_lc = fill(T(-1),length(c_to_a))
  for (c,a) in enumerate(c_to_a)
    b = c_to_b[c]
    for (lc,p) in enumerate(datarange(b_to_la_to_a,b))
      if a == b_to_la_to_a.data[p]
        c_to_lc[c] = T(lc)
        break
      end
    end
  end
  return c_to_lc
end

function find_local_index(
  ::Type{T}, c_to_la_to_a :: AbstractVector{<:AbstractVector}, c_to_b :: AbstractVector, b_to_la_to_a :: Table
) where T
  c1 = array_cache(c_to_la_to_a)
  c2 = array_cache(b_to_la_to_a)
  ptrs = generate_ptrs(c_to_la_to_a)
  data = fill(T(-1),ptrs[end]-1)
  for c in eachindex(c_to_la_to_a)
    b = c_to_b[c]
    lin_to_a = getindex!(c1,c_to_la_to_a,c)
    lout_to_a = getindex!(c2,b_to_la_to_a,b)

    pin, pout = sortperm(lin_to_a), sortperm(lout_to_a)
    nin, nout = length(lin_to_a), length(lout_to_a)
    kin, kout = 1, 1
    while kin <= nin && kout <= nout
      ain = lin_to_a[pin[kin]]
      aout = lout_to_a[pout[kout]]
      if ain == aout
        k = ptrs[c] + pin[kin] - 1
        data[k] = T(pout[kout])
        kin += 1
        kout += 1
      elseif ain < aout
        kin += 1
      else
        kout += 1
      end
    end
  end
  return Table(data,ptrs)
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

function flatten_partition!(b_to_a,a_to_bs::Table)
  for a in eachindex(a_to_bs)
    for p in datarange(a_to_bs,a)
      b = a_to_bs.data[p]
      b_to_a[b] = a
    end
  end
end

"""
    find_local_nbor_index(a_to_b, a_to_lb_to_b) -> a_to_lb
"""
function find_local_nbor_index(a_to_b, a_to_lb_to_b::Table)
  a_to_lb = Vector{GridapLocalInt}(undef,length(a_to_b))
  for (a,b) in enumerate(a_to_b)
    for (lb,p) in enumerate(datarange(a_to_lb_to_b,a))
      if b == a_to_lb_to_b.data[p]
        a_to_lb[a] = GridapLocalInt(lb)
        break
      end
    end
  end
  return a_to_lb
end

"""
    find_local_nbor_index(a_to_b, a_to_c, c_to_lb_to_b) -> a_to_lb
"""
function find_local_nbor_index(a_to_b, a_to_c, c_to_lb_to_b::Table)
  a_to_lb = Vector{GridapLocalInt}(undef,length(a_to_b))
  for (a,b) in enumerate(a_to_b)
    c = a_to_c[a]
    for (lb,p) in enumerate(datarange(c_to_lb_to_b,c))
      if b == c_to_lb_to_b.data[p]
        a_to_lb[a] = GridapLocalInt(lb)
        break
      end
    end
  end
  return a_to_lb
end

"""
    merge_entries(a_to_lb_to_b, c_to_la_to_a) -> c_to_lb_to_b

Merge the entries of `a_to_lb_to_b`, grouping them by `c_to_la_to_a`. Returns 
the merged table `c_to_lb_to_b`.

Accepts the following keyword arguments:

- `acc`: Accumulator for the entries of `a_to_lb_to_b`. Default to a `Set`, ensuring 
         that the resulting entries are unique.
- `post`: Postprocessing function to apply to the accumulator before storing the resulting entries.
          Defaults to the identity, but can be used to perform local sorts or filters, for example.
"""
function merge_entries(
  a_to_lb_to_b::AbstractVector{<:AbstractVector{T}},
  c_to_la_to_a::AbstractVector{<:AbstractVector{Ti}}; 
  acc  = Set{T}(),
  post = identity
) where {T,Ti<:Integer}
  c1 = array_cache(a_to_lb_to_b)
  c2 = array_cache(c_to_la_to_a)

  n_c = length(c_to_la_to_a)
  ptrs = zeros(Int32,n_c+1)
  for c in 1:n_c
    as = getindex!(c2,c_to_la_to_a,c)
    for a in as
      bs = getindex!(c1,a_to_lb_to_b,a)
      !isempty(bs) && push!(acc, bs...)
    end
    ptrs[c+1] += length(post(acc))
    empty!(acc)
  end
  length_to_ptrs!(ptrs)

  data = zeros(T,ptrs[end]-1)
  for c in 1:n_c
    as = getindex!(c2,c_to_la_to_a,c)
    for a in as
      bs = getindex!(c1,a_to_lb_to_b,a)
      !isempty(bs) && push!(acc, bs...)
    end
    data[ptrs[c]:ptrs[c+1]-1] = post(collect(acc))
    empty!(acc)
  end

  c_to_lb_to_b = Table(data,ptrs)
  return c_to_lb_to_b
end

"""
    block_identity_array(::Type{T},ptrs) where T

Given a vector of pointers of length `n+1`, returns a vector of length `ptrs[end]-1` 
where the entries are the index of the block to which each entry belongs.

# Example

```julia

julia> block_identity_array([1,3,7])

6-element Vector{Int64}:
 1
 1
 2
 2
 2
 2
```
"""
function block_identity_array(::Type{T},ptrs) where T
  n = length(ptrs)-1
  a = Vector{T}(undef,ptrs[end]-1)
  for i in 1:n
    a[ptrs[i]:ptrs[i+1]-1] .= i
  end
  return a
end

block_identity_array(ptrs) = block_identity_array(Int,ptrs)

"""
    local_identity_array(::Type{T}, ptrs) where T

Given a vector of pointers of length `n+1`, returns a vector of length `ptrs[end]-1`
where the entries are the local index of the entry within the block it belongs to.

# Example

```julia

julia> local_identity_array([1,3,7])

6-element Vector{Int64}:
 1
 2
 1
 2
 3
 4
```
"""
function local_identity_array(::Type{T}, ptrs) where T
  n = length(ptrs)-1
  a = Vector{T}(undef,ptrs[end]-1)
  for i in 1:n
    ni = ptrs[i+1]-ptrs[i]
    a[ptrs[i]:ptrs[i+1]-1] .= 1:ni
  end
  return a
end

local_identity_array(ptrs) = local_identity_array(Int,ptrs)

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
