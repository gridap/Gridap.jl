
struct VectorOfBlockArrayCoo{T,N,B,X,Z} <: AbstractVector{T}
  blocks::B
  blockids::Vector{NTuple{N,Int}}
  axes::X
  ptrs::Array{Int,N}
  zero_blocks::Z
  function VectorOfBlockArrayCoo(
    blocks::Tuple,
    blockids::Vector{NTuple{N,Int}},
    axes::AbstractArray{<:NTuple{N}},
    ptrs::Array{Int,N} = _compute_ptrs(blockids,testitem(axes)),
    zero_blocks::Tuple=_compute_zero_blocks_array(blocks,ptrs,axes)) where N

    msg = "Trying to build a VectorOfBlockArrayCoo with repeated blocks"
    @assert _no_has_repeaded_blocks(blockids,ptrs) msg

    if length(blocks) != 0
      @assert all(map(length,blocks) .==  length(first(blocks)) )
      blocks_i = collect(map(testitem,blocks))
      zero_blocks_i = collect(eltype(blocks_i),map(testitem,zero_blocks))
    else
      @assert length(zero_blocks) !=0
      zero_blocks_i = collect(map(testitem,zero_blocks))
      blocks_i = collect(eltype(zero_blocks_i),map(testitem,blocks))
    end

    B = typeof(blocks)
    X = typeof(axes)
    axes_i = testitem(axes)
    t = BlockArrayCoo(blocks_i,blockids,axes_i,ptrs,zero_blocks_i)
    T = typeof(t)
    Z = typeof(zero_blocks)
    new{T,N,B,X,Z}(blocks,blockids,axes,ptrs,zero_blocks)
  end
end

function VectorOfBlockArrayCoo(
  blocks::Tuple,
  blockids::Vector{NTuple{N,Int}},
  axes::AbstractArray{<:NTuple{N}},
  zero_blocks::Tuple) where N

  ptrs = _compute_ptrs(blockids,testitem(axes))
  VectorOfBlockArrayCoo(blocks,blockids,axes,ptrs,zero_blocks)
end

const VectorOfBlockVectorCoo = VectorOfBlockArrayCoo{T,1} where T
const VectorOfBlockMatrixCoo = VectorOfBlockArrayCoo{T,2} where T

function _no_has_repeaded_blocks(blockids::Vector{NTuple{N,Int}},ptrs) where N
  maxblocks = size(ptrs)
  touched = zeros(Int,maxblocks)
  for c in blockids
    touched[c...] += 1
  end
  all( touched .<= 1 )
end

function _compute_zero_blocks_array(blocks,ptrs,axs)
  A = eltype(first(blocks))
  cis = CartesianIndices(ptrs)
  zero_blocks = []
  for ci in cis
    p = ptrs[ci]
    if p<0
      block = _zeros_like(A,Tuple(ci),axs)
      push!(zero_blocks,block)
    end
  end
  Tuple(zero_blocks)
end

function zeros_like(a::VectorOfBlockArrayCoo{T,N} where T) where N
  ptrs = copy(a.ptrs)
  for k in 1:length(ptrs)
    ptrs[k] = -k
  end
  zero_blocks = _compute_zero_blocks_array(a.blocks,ptrs,a.axes)
  blocks = ()
  blockids = NTuple{N,Int}[]
  VectorOfBlockArrayCoo(blocks,blockids,a.axes,ptrs,zero_blocks)
end

function _zeros_like(::Type{A},I,axs) where A
  block = apply(axs) do a
    laxs = map( local_range, a, I)
    zeros_like(A,laxs)
  end
  block
end

function _zeros_like(::Type{<:BlockArrayCoo{T,N,A} where T},I,axs) where {N,A}
  axsI = apply(axs) do a
    map(local_range, a, I)
  end
  blocks = ()
  blockids = NTuple{N,Int}[]
  ptrs = zeros(Int,map(x->blocksize(x,1),testitem(axsI)))
  for k in 1:length(ptrs)
    ptrs[k] = -k
  end
  zero_blocks = []
  cis = CartesianIndices(ptrs)
  for ci in cis
    block = _zeros_like(A,Tuple(ci),axsI)
    push!(zero_blocks,block)
  end
  VectorOfBlockArrayCoo(blocks,blockids,axsI,ptrs,Tuple(zero_blocks))
end

BlockArrays.blocksize(a::VectorOfBlockArrayCoo) = blocksize(testitem(a))

BlockArrays.blocksize(a::VectorOfBlockArrayCoo,i::Integer) = blocksize(testitem(a),i)

function BlockArrays.eachblock(a::VectorOfBlockArrayCoo)
  cis = CartesianIndices(blocksize(a))
  blocks = map(ci->Block(Tuple(ci)),cis)
  ( a[block] for block in blocks  )
end

Base.IndexStyle(::Type{<:VectorOfBlockArrayCoo}) = IndexLinear()

function Base.size(a::VectorOfBlockArrayCoo) 
  if length(a.blocks) != 0
    (length(first(a.blocks)),)
  else
    (length(first(a.zero_blocks)),)
  end
end

function array_cache(a::VectorOfBlockArrayCoo)

  if length(a.blocks) != 0
    blocks_i = collect(map(testitem,a.blocks))
    zero_blocks_i = collect(eltype(blocks_i),map(testitem,a.zero_blocks))
  else
    zero_blocks_i = collect(map(testitem,a.zero_blocks))
    blocks_i = collect(eltype(zero_blocks_i),map(testitem,a.blocks))
  end

  ca = array_cache(a.axes)
  cb = array_caches(a.blocks...)
  cz = array_caches(a.zero_blocks...)
  (blocks_i,zero_blocks_i,ca,cb,cz)
end

@inline function getindex!(cache,a::VectorOfBlockArrayCoo,i::Integer)
  blocks_i, zero_blocks_i, ca, cb, cz = cache
  axes_i = getindex!(ca,a.axes,i)
  blocks_i .= getitems!(cb,a.blocks,i)
  zero_blocks_i .= getitems!(cz,a.zero_blocks,i)
  BlockArrayCoo(blocks_i,a.blockids,axes_i,a.ptrs,zero_blocks_i)
end

function Base.getindex(a::VectorOfBlockArrayCoo,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function Base.getindex(a::VectorOfBlockArrayCoo,b::Block)
  i = convert(Tuple,b)
  p = a.ptrs[i...]
  if p>0
    a.blocks[p]
  else
    a.zero_blocks[-p]
  end
end

function Base.getindex(a::VectorOfBlockArrayCoo,b::Block{1}...)
  a[Block(map(Int,b)...)]
end

function is_zero_block(a::VectorOfBlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  @assert p != 0
  p < 0
end

function apply(::typeof(transpose),a::VectorOfBlockMatrixCoo)
  blocks = [ apply(transpose,block) for block in a.blocks ]
  zero_blocks = [ apply(transpose,block) for block in a.zero_blocks ]
  blockids = [ (j,i) for (i,j) in a.blockids ]
  axs = apply( ax->(ax[2],ax[1]), a.axes)
  ptrs = collect(Transpose(a.ptrs))
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs,ptrs,Tuple(zero_blocks))
end

