module BlockArraysCOOTests


using Gridap.Helpers
import LinearAlgebra: mul!
using Gridap.Arrays

"""
"""
struct BlockArrayCOO{T,N,A<:AbstractArray{T,N}} <: GridapType
  blocks::Vector{A}
  coordinates::Vector{NTuple{N,Int}}
  ptrs::Array{Int,N}
  block_size::NTuple{N,Int}

  function BlockArrayCOO(
    blocks::Vector{A},
    coordinates::Vector{NTuple{N,Int}}) where {T,N,A<:AbstractArray{T,N}}

    @assert length(blocks) == length(coordinates)
    msg = "Trying to build a BlockArrayCOO with repeated blocks"
    @assert _no_has_repeaded_blocks(coordinates) msg
    ptrs = _prepare_ptrs(coordinates)
    block_size = _get_block_size(coordinates)
    new{T,N,A}(blocks,coordinates,ptrs,block_size)
  end
end

function _no_has_repeaded_blocks(coordinates::Vector{NTuple{N,Int}}) where N
  maxblocks = _get_block_size(coordinates)
  touched = zeros(Int,maxblocks)
  for c in coordinates
    touched[c...] += 1
  end
  all( touched .<= 1 )
end

function _prepare_ptrs(coordinates)
  s = _get_block_size(coordinates)
  ptrs = zeros(Int,s)
  for (i,c) in enumerate(coordinates)
    ptrs[c...] = i
  end
  ptrs
end

function _get_block_size(coordinates::Vector{NTuple{N,Int}}) where N
  m = zeros(Int,N)
  for c in coordinates
    for n in 1:N
      m[n] = max(m[n],c[n])
    end
  end
  NTuple{N,Int}(m)
end

function get_block_size(a::BlockArrayCOO)
  a.block_size
end

function num_blocks(a::BlockArrayCOO)
  s = get_block_size(a)
  prod(s)
end

function num_stored_blocks(a::BlockArrayCOO)
  length(a.coordinates)
end

function has_all_blocks(a::BlockArrayCOO{T,N}) where {T,N}
  num_blocks(a) == num_stored_blocks(a)
end

"""
"""
function add_to_array!(a::AbstractArray{Ta,N},b::AbstractArray{Tb,N},combine=+) where {Ta,Tb,N}
  @assert size(a) == size(b) "Arrays sizes mismatch"
  @inbounds for i in eachindex(a)
    a[i] = combine(a[i],b[i])
  end
end

"""
"""
function add_to_array!(a::AbstractArray,b::Number,combine=+)
  @inbounds for i in eachindex(a)
    a[i] = combine(a[i],b)
  end
end

function Base.:*(a::BlockArrayCOO,b::Number)
  blocks = [ block*b for block in a.blocks ]
  coordinates = a.coordinates
  BlockArrayCOO(blocks,coordinates)
end

function Base.:*(b::Number,a::BlockArrayCOO)
  blocks = [ b*block for block in a.blocks ]
  coordinates = a.coordinates
  BlockArrayCOO(blocks,coordinates)
end

function Base.show(io::IO,::MIME"text/plain",a::BlockArrayCOO)
  println(io,"BlockArrayCOO object:")
  cis = CartesianIndices(a.ptrs)
  for ci in cis
    p = a.ptrs[ci]
    if p == 0
      println(io,"Block $(Tuple(ci)) -> Empty")
    else
      println(io,"Block $(Tuple(ci)) -> $(a.blocks[p])")
    end
  end
end

function add_to_array!(a::BlockArrayCOO{Ta,N},b::BlockArrayCOO{Tb,N},combine=+) where {Ta,Tb,N}
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    bk = b.blocks[k]
    add_to_array!(ak,bk,combine)
  end
end

function add_to_array!(a::BlockArrayCOO,b::Number,combine=+)
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    add_to_array!(ak,b,combine)
  end
end

function Base.:*(a::BlockArrayCOO{Ta,2},b::BlockArrayCOO{Tb,1}) where {Ta,Tb}
  @assert num_stored_blocks(a) != 0
  @notimplementedif ! has_all_blocks(b)

  function fun(i,a,b)
    ai = a.blocks[i]
    ci, cj = a.coordinates[i]
    p = b.ptrs[cj]
    bi = b.blocks[p]
    ai*bi
  end

  blocks = [ fun(i,a,b) for i in 1:length(a.blocks)]
  coordinates = [ (c[1],)  for c in a.coordinates]

  data = _merge_repeated_blocks(blocks,coordinates)
  BlockArrayCOO(data...)
end

function _merge_repeated_blocks(blocks,coordinates::Vector{NTuple{N,Int}}) where N
  @assert length(blocks) == length(coordinates)
  s = _get_block_size(coordinates)
  ptrs = zeros(Int,s)
  A = eltype(blocks)
  _blocks = A[]
  _coords = NTuple{N,Int}[]
  p = 1
  for b in 1:length(blocks)
    c = coordinates[b]
    block = blocks[b]
    p = ptrs[c...]
    if p == 0
      push!(_blocks,block)
      push!(_coords,c)
      p += 1
      ptrs[c...] = p
    else
      add_to_array!(_blocks[p],block)
    end
  end
  (_blocks,_coords)
end

function mul!(c::BlockArrayCOO{Tc,1},a::BlockArrayCOO{Ta,2},b::BlockArrayCOO{Tb,1}) where {Tc,Ta,Tb}
  for ci in c.blocks
    fill!(ci,zero(Tc))
  end
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    ci, cj = a.coordinates[k]
    p = b.ptrs[cj]
    bk = b.blocks[p]
    q = c.ptrs[ci]
    ck = c.blocks[q]
    mul!(ck,ak,bk,1,1)
  end
end

function CachedBlockArrayCOO(a::BlockArrayCOO)
  blocks = [ CachedArray(b) for b in a.blocks ]
  coordinates = a.coordinates
  BlockArrayCOO(blocks,coordinates)
end

function _resize_for_mul!(
  c::BlockArrayCOO{Tc,1},a::BlockArrayCOO{Ta,2},b::BlockArrayCOO{Tb,1}) where {Tc,Ta,Tb}
  
  for k in 1:length(a.blocks)
    ak = a.blocks[k]
    ci, cj = a.coordinates[k]
    q = c.ptrs[ci]
    ck = c.blocks[q]
    setsize!(ck,(size(ak,1),))
  end

end

function _move_cached_arrays!(r::BlockArrayCOO,c::BlockArrayCOO)
  for  k in 1:length(c.blocks)
    ck = c.blocks[k]
    r.blocks[k] = ck.array
  end
end

function Base.getindex(a::BlockArrayCOO{T,N},I::Vararg{Int,N}) where {T,N}
  p = a.ptrs[I...]
  @assert p > 0 "You are attempting to access a block that is not stored"
  a.blocks[p]
end

using Test

A11 = ones(Int,2,3)
A21 = 2*ones(Int,4,3)
A12 = 3*ones(Int,2,5)
blocks = [A11,A21,A12]
coordinates = [(1,1),(2,1),(1,2)]
a = BlockArrayCOO(blocks,coordinates)

@test a.blocks == blocks
@test a.coordinates == coordinates
@test a.ptrs == [1 3; 2 0]
@test num_blocks(a) == 4
@test num_stored_blocks(a) == 3
@test get_block_size(a) == (2,2)
@test has_all_blocks(a) == false
@test a[1,1] == A11
@test a[2,1] == A21

B1 = 10*ones(Int,3)
B2 = 20*ones(Int,5)
blocks = [B1,B2]
coordinates = [(1,),(2,)]
b = BlockArrayCOO(blocks,coordinates)

c = a*b

@test c.blocks[1] == A11*B1 + A12*B2
@test c.blocks[2] == A21*B1
@test c.coordinates == [(1,),(2,)]

mul!(c,a,b)

@test c.blocks[1] == A11*B1 + A12*B2
@test c.blocks[2] == A21*B1
@test c.coordinates == [(1,),(2,)]

c = a*b
r = CachedBlockArrayCOO(c)

_resize_for_mul!(r,a,b)
_move_cached_arrays!(c,r)
mul!(c,10*a,b)

_resize_for_mul!(r,a,b)
_move_cached_arrays!(c,r)
mul!(c,a,b)

A11 = ones(Int,7,3)
A21 = 2*ones(Int,4,3)
A12 = 3*ones(Int,7,5)
blocks = [A11,A21,A12]
coordinates = [(1,1),(2,1),(1,2)]
_a = BlockArrayCOO(blocks,coordinates)

_resize_for_mul!(r,_a,b)
_move_cached_arrays!(c,r)
mul!(c,_a,b)

_resize_for_mul!(r,a,b)
_move_cached_arrays!(c,r)
mul!(c,a,b)

add_to_array!(c,c)
add_to_array!(c,10)

@show c
display(a)

end # module
