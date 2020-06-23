
struct MultiFieldCellArray{T,N,B<:Tuple} <: AbstractVector{MultiFieldArray{T,N,Array{T,N}}}
  blocks::B
  block_ids::Vector{NTuple{N,Int}}
  function MultiFieldCellArray(_blocks::Tuple,_block_ids::Vector{NTuple{N,Int}}) where N
    blocks, block_ids = _merge_repeated_blocks_mca(_blocks,_block_ids)
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

function _merge_repeated_blocks_mca(blocks,coordinates::Vector{NTuple{N,Int}}) where N
  @assert length(blocks) == length(coordinates)
  s = _get_block_size(coordinates)
  ptrs = zeros(Int,s)
  _blocks = []
  _coords = NTuple{N,Int}[]
  q = 1
  for b in 1:length(blocks)
    c = coordinates[b]
    block = blocks[b]
    p = ptrs[c...]
    if p == 0
      push!(_blocks,block)
      push!(_coords,c)
      ptrs[c...] = q
      q += 1
    else
      _blocks[p] = apply(elem(+),_blocks[p],block)
    end
  end
  ( tuple(_blocks...), _coords)
end

function array_cache(a::MultiFieldCellArray{T,N}) where {T,N}
  coordinates = a.block_ids
  nblocks = length(coordinates)
  blocks = Vector{Array{T,N}}(undef,nblocks)
  b = MultiFieldArray(blocks,coordinates)
  caches = array_caches(a.blocks...)
  (b,caches)
end

function getindex!(cache,a::MultiFieldCellArray,i::Integer)
  b, caches = cache
  bis = getitems!(caches,a.blocks,i)
  for (k,bk) in enumerate(bis)
    b.blocks[k] = bk
  end
  b
end

function Base.getindex(a::MultiFieldCellArray,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function Base.size(a::MultiFieldCellArray) 
  bi, = a.blocks
  size(bi)
end

function reindex(a::MultiFieldCellArray,b::AbstractArray)
  f = (ai) -> reindex(ai,b)
  blocks = map(f,a.blocks)
  MultiFieldCellArray(blocks,a.block_ids)
end

function reindex(a::MultiFieldCellArray,b::IdentityVector)
  a
end

struct BlockTracker{N} <: GridapType
  blocks::Tuple
  block_ids::Vector{NTuple{N,Int}}
end

function operate(::typeof(+),a::BlockTracker)
  a
end

function operate(::typeof(-),a::BlockTracker)
  new_blocks = map(-,a.blocks)
  BlockTracker(new_blocks,a.block_ids)
end

function operate(::typeof(+),a::BlockTracker,b::BlockTracker)
  new_blocks = (a.blocks...,b.blocks...)
  new_block_ids = vcat(a.block_ids,b.block_ids)
  BlockTracker(new_blocks,new_block_ids)
end

function operate(::typeof(-),a::BlockTracker,b::BlockTracker)
  new_blocks = (a.blocks...,(-b).blocks...)
  new_block_ids = vcat(a.block_ids,b.block_ids)
  BlockTracker(new_blocks,new_block_ids)
end

function operate(op,a::BlockTracker,b::BlockTracker)
  msg = "Operation $op not yet implemented in this context"
  @notimplementedif !( op  in (*,inner,dot) ) msg
  new_blocks = []
  new_block_ids = NTuple{2,Int}[]
  for i in 1:length(a.blocks)
    ai = a.blocks[i]
    @notimplementedif ! isa(ai,CellBasisWithFieldID) msg
    ai_id, = a.block_ids[i]
    for j in 1:length(b.blocks)
      bj = b.blocks[j]
      @notimplementedif ! isa(bj,CellBasisWithFieldID) msg
      bj_id, = b.block_ids[j]
      push!(new_blocks, operate(op,ai.cell_basis,bj.cell_basis))
      push!(new_block_ids, (ai_id,bj_id))
    end
  end
  BlockTracker(Tuple(new_blocks),new_block_ids)
end

function operate(op,a,b::BlockTracker)
  _operate_bt_b(op,a,b)
end

function operate(op,a::CellField,b::BlockTracker)
  _operate_bt_b(op,a,b)
end

function operate(op,b::BlockTracker,a)
  _operate_bt_a(op,b,a)
end

function operate(op,b::BlockTracker,a::CellField)
  _operate_bt_a(op,b,a)
end

function _operate_bt_b(op,a,b)
  msg = "Operation $op not yet implemented in this context"
  @notimplementedif !( op  in (*,inner) ) msg
  new_blocks = map(bi->operate(op,a,bi),b.blocks)
  BlockTracker(new_blocks,b.block_ids)
end

function _operate_bt_a(op,b,a)
  msg = "Operation $op not yet implemented in this context"
  @notimplementedif !( op  in (*,inner) ) msg
  new_blocks = map(bi->operate(op,a,bi),b.blocks)
  BlockTracker(new_blocks,b.block_ids)
end

