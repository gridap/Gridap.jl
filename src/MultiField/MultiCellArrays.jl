
struct MultiCellArray{T,N,B<:Tuple} <: AbstractVector{BlockArrayCOO{T,N,Array{T,N}}}
  blocks::B
  block_ids::Vector{NTuple{N,Int}}
  function MultiCellArray(blocks::Tuple,block_ids::Vector{NTuple{N,Int}}) where N
    @assert length(blocks) > 0
    @assert length(blocks) == length(block_ids)
    bi, = blocks
    @assert isa(bi,AbstractArray)
    @assert all( ( size(b) == size(bi) for b in blocks ) )
    @assert all( ( eltype(b) == eltype(bi) for b in blocks ) )
    A = eltype(bi)
    @assert A <: Array
    @assert ndims(A) == N
    T = eltype(A)
    B = typeof(blocks)
    new{T,N,B}(blocks,block_ids)
  end
end

function array_cache(a::MultiCellArray{T,N}) where {T,N}
  coordinates = a.block_ids
  nblocks = length(coordinates)
  blocks = Vector{Array{T,N}}(undef,nblocks)
  b = BlockArrayCOO(blocks,coordinates)
  caches = array_caches(a.blocks...)
  (b,caches)
end

function getindex!(cache,a::MultiCellArray,i::Integer)
  b, caches = cache
  bis = getitems!(caches,a.blocks,i)
  for (k,bk) in enumerate(bis)
    b.blocks[k] = bk
  end
  b
end

function Base.getindex(a::MultiCellArray,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function Base.size(a::MultiCellArray) 
  bi, = a.blocks
  size(bi)
end

function reindex(a::MultiCellArray,b::AbstractArray)
  f = (ai) -> reindex(ai,b)
  blocks = map(f,a.blocks)
  MultiCellArray(blocks,a.block_ids)
end

function reindex(a::MultiCellArray,b::IdentityVector)
  a
end

struct BlockTracker{N} <: GridapType
  blocks::Tuple
  block_ids::Vector{NTuple{N,Int}}
end

function Base.:+(a::BlockTracker)
  a
end

function Base.:-(a::BlockTracker)
  new_blocks = map(-,a.blocks)
  BlockTracker(new_blocks,a.block_ids)
end

function Base.:+(a::BlockTracker,b::BlockTracker)
  new_blocks = (a.blocks...,b.blocks...)
  new_block_ids = vcat(a.block_ids,b.block_ids)
  BlockTracker(new_blocks,new_block_ids)
end

function Base.:-(a::BlockTracker,b::BlockTracker)
  new_blocks = (a.blocks...,(-b).blocks...)
  new_block_ids = vcat(a.block_ids,b.block_ids)
  BlockTracker(new_blocks,new_block_ids)
end


