
struct BlockArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  blocks::Vector{A}
  blockids::Vector{NTuple{N,Int}}
  axes::X
  ptrs::Array{Int,N}
  zero_blocks::Vector{A}
  function BlockArrayCoo(
    blocks::Vector{A},
    blockids::Vector{NTuple{N,Int}},
    axes::NTuple{N},
    ptrs::Array{Int,N}=_compute_ptrs(blockids,axes),
    zero_blocks::Vector{A}=_compute_zero_blocks(A,ptrs,axes)) where {T,N,A<:AbstractArray{T,N}}

    X = typeof(axes)
    new{T,N,A,X}(blocks,blockids,axes,ptrs,zero_blocks)
  end
end

function _compute_ptrs(blockids,axes)
  s = map(i->first(blocksize(i)),axes)
  ptrs = zeros(Int,s)
  for (i,c) in enumerate(blockids)
    ptrs[c...] = i
  end
  m = 1
  for (k,p) in enumerate(ptrs)
    if p ==0
      ptrs[k] = -m
      m += 1
    end
  end
  ptrs
end

function _compute_zero_blocks(::Type{A},ptrs,axes) where A
  cis = CartesianIndices(ptrs)
  zero_blocks = A[]
  for ci in cis
    p = ptrs[ci]
    if p<0
      i = Tuple(ci)
      s = map((k,l)->length(k[Block(l)]),axes,i)
      block = A(undef,s)
      fill!(block,zero(eltype(A)))
      push!(zero_blocks,block)
    end
  end
  zero_blocks
end

function BlockArrays.getblock(a::BlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  if p>0
    a.blocks[p]
  else
    a.zero_blocks[-p]
  end
end

@inline is_nonzero_block(a,b...) = ! is_zero_block(a,b...)

function is_zero_block(a,b::Block)
  i = convert(Tuple,b)
  is_zero_block(a,i...)
end

function is_zero_block(a,b::Block...)
  i = map(Int,b)
  is_zero_block(a,i...)
end

function is_zero_block(a,b::CartesianIndex)
  i = Tuple(b)
  is_zero_block(a,i...)
end

function is_zero_block(a::BlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  @assert p != 0
  p < 0
end

function BlockArrays.getblock!(block,a::BlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  if p>0
    block .= a.blocks[p]
  else
    block .= a.zero_blocks[-p]
  end
end

function BlockArrays.eachblock(a::BlockArrayCoo)
  cis = CartesianIndices(blocksize(a))
  blocks = map(ci->Block(Tuple(ci)),cis)
  ( a[block] for block in blocks  )
end

function enumerateblocks(a)
  cis = CartesianIndices(blocksize(a))
  blocks = map(ci->Block(Tuple(ci)),cis)
  zip(blocks,eachblock(a))
end

Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))

Base.axes(a::BlockArrayCoo) = a.axes

function Base.getindex(a::BlockArrayCoo,i::Integer...)
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

