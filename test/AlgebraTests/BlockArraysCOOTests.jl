module BlockArraysCOOTests


using Gridap.Helpers
import LinearAlgebra: mul!

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

"""
"""
function get_block_size(a::BlockArrayCOO)
  a.block_size
end

"""
"""
function num_blocks(a::BlockArrayCOO)
  s = get_block_size(a)
  prod(s)
end

"""
"""
function num_stored_blocks(a::BlockArrayCOO)
  length(a.coordinates)
end

"""
"""
function has_all_blocks(a::BlockArrayCOO{T,N}) where {T,N}
  num_blocks(a) == num_stored_blocks(a)
end

#function Base.:*(a::BlockArrayCOO{Ta,2},b::BlockArrayCOO{Tb,1}) where {Ta,Tb}
#
#
#end

"""
"""
function add_to_array!(a::AbstractArray{Ta,N},b::AbstractArray{Tb,N}) where {Ta,Tb,N}
  @assert size(a) == size(b) "Arrays sizes mismatch"
  @inbounds for i in eachindex(a)
    a[i] += b[i]
  end
end

"""
"""
function add_to_array!(a::AbstractArray,b::Number)
  @inbounds for i in eachindex(a)
    a[i] += b
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

end # module
