module CellMaps

using Test
using Gridap
using Gridap.Helpers
using Gridap.CachedArrays
using Gridap.CellValues: _test_iter_cell_value
using Gridap.CellValues: _test_index_cell_value
using Gridap.CellValues: _eq

export CellMap
export IterCellMap
export IndexCellMap
import Gridap: evaluate
export test_iter_cell_map
export test_index_cell_map
export test_index_cell_map_with_index_arg
export test_iter_cell_map_with_index_result

import Base: iterate
import Base: length
import Base: size
import Base: getindex

"""
Abstract object that traverses a set of cells and at every cell returns a
`Map{S,M,T,N}`
"""
const IterCellMap{S,M,T,N,R<:Map{S,M,T,N}} = IterCellValue{R}

"""
Abstract array indexed by cells that returns a `Map{S,M,T,N}`
"""
const IndexCellMap{S,M,T,N,C,R<:Map{S,M,T,N}} = IndexCellValue{R,C}

"""
Abstract object that for a given cell index returns a `Map{S,M,T,N}`
"""
const CellMap{S,M,T,N} = Union{IterCellMap{S,M,T,N},IndexCellMap{S,M,T,N}}

# evaluation

"""
Return the cellwise maps of a `CellMap` on a cellwise set of points
"""
function evaluate(cm::CellMap{S,M,T,N},ca::CellArray{<:S,M}) where {S,M,T,N}
  IterCellMapValue(cm,ca)
end

function evaluate(cm::IndexCellMap{S,M,T,N},ca::IndexCellArray{<:S,M}) where {S,M,T,N}
  IndexCellMapValue(cm,ca)
end

struct IterCellMapValue{T,N,V,A} <: IterCellArray{T,N,CachedArray{T,N,Array{T,N}}}
  cm::V
  ca::A
end

function IterCellMapValue(cm::CellMap{S,M,T,N},ca::CellArray{<:S,M}) where {S,M,T,N}
  V = typeof(cm)
  A = typeof(ca)
  IterCellMapValue{T,N,V,A}(cm,ca)
end

length(v::IterCellMapValue) = length(v.ca)

@inline function iterate(v::IterCellMapValue{T,N}) where {T,N}
  cmnext = iterate(v.cm)
  canext = iterate(v.ca)
  r = CachedArray(T,N)
  _iterate(cmnext,canext,r)
end

@inline function iterate(v::IterCellMapValue,state)
  cmstate, castate, r = state
  cmnext = iterate(v.cm,cmstate)
  canext = iterate(v.ca,castate)
  _iterate(cmnext,canext,r)
end

@inline function _iterate(cmnext,canext,r)
  if cmnext === nothing; return nothing end
  if canext === nothing; return nothing end
  m, cmstate = cmnext
  a, castate = canext
  s = return_size(m,size(a))
  setsize!(r,s)
  evaluate!(m,a,r)
  state = (cmstate,castate,r)
  (r, state)
end

struct IndexCellMapValue{T,N,V,A} <: IndexCellArray{T,N,CachedArray{T,N,Array{T,N}},1}
  cm::V
  ca::A
  r::CachedArray{T,N,Array{T,N}}
end

function IndexCellMapValue(cm::IndexCellMap{S,M,T,N},ca::IndexCellArray{<:S,M}) where {S,M,T,N}
  V = typeof(cm)
  A = typeof(ca)
  r = CachedArray(T,N)
  IndexCellMapValue{T,N,V,A}(cm,ca,r)
end

length(v::IndexCellMapValue) = length(v.ca)

size(v::IndexCellMapValue) = (length(v),)

function getindex(v::IndexCellMapValue,i::Integer)
  m = v.cm[i]
  a = v.ca[i]
  s = return_size(m,size(a))
  setsize!(v.r,s)
  evaluate!(m,a,v.r)
  v.r
end

# Testers

function test_iter_cell_map(
  m::CellMap{S,M,T,N},
  a::CellArray{<:S,M},
  b::AbstractArray{<:AbstractArray{T,N}}) where {S,M,T,N}

  @test length(m) == length(b)
  @test length(m) == length(a)

  c = evaluate(m,a)
  _test_iter_cell_value(c,b)

  for (mi,ai,bi) in zip(m,a,b)
    @assert isa(mi,Map{S,M,T,N})
    @assert _eq(evaluate(mi,ai),bi)
    @assert typeof(mi) == eltype(m)
  end

end

function test_index_cell_map(
  m::IndexCellMap{S,M,T,N},
  a::CellArray{<:S,M},
  b::AbstractArray{<:AbstractArray{T,N}}) where {S,M,T,N}

  @test size(m) == size(b)

  @test IndexStyle(m) == IndexStyle(b)

  test_iter_cell_map(m,a,b)

  for (i,ai) in enumerate(a)
    mi = m[i]
    bi = b[i]
    @assert isa(mi,Map{S,M,T,N})
    @assert _eq(evaluate(mi,ai),bi)
    @assert typeof(mi) == eltype(m)
  end

end

function test_index_cell_map_with_index_arg(
  m::IndexCellMap{S,M,T,N},
  a::IndexCellArray{<:S,M},
  b::AbstractArray{<:AbstractArray{T,N}}) where {S,M,T,N}

  test_index_cell_map(m,a,b)

  c = evaluate(m,a)
  @test isa(c,IndexCellArray)

  @test length(c) == length(m)
  @test length(c) == length(a)

  for i in 1:length(c)
    mi = m[i]
    bi = b[i]
    ci = c[i]
    @assert _eq(bi,ci)
  end

  for i in 1:length(m)
    mi = m[i]
    ai = a[i]
    bi = b[i]
    @assert isa(mi,Map{S,M,T,N})
    @assert _eq(evaluate(mi,ai),bi)
    @assert typeof(mi) == eltype(m)
  end

end

function test_iter_cell_map_with_index_result(
  m::CellMap{S,M,T,N},
  a::CellArray{<:S,M},
  b::AbstractArray{<:AbstractArray{T,N}}) where {S,M,T,N}

  test_iter_cell_map(m,a,b)

  c = evaluate(m,a)

  _test_index_cell_value(c,b)

end

end # module CellMaps
