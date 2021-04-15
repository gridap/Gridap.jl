
# Restricted for a single non zero block
struct BlockFieldArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  axes::X
  blockids::Vector{NTuple{N,Int}}
  block::A
  function BlockFieldArrayCoo(
    _axes::NTuple{N},
    blockids::Vector{NTuple{N,Int}},
    block::A) where {T,N,A<:AbstractArray{T,N}}

    @assert length(blockids) == 1
    #I = first(blockids)
    #@check blocks_equal(axes(block),map(local_range,_axes,I)) "The given block and axes are incompatible."

    X = typeof(_axes)
    new{T,N,A,X}(_axes,blockids,block)
  end
end

struct BlockFieldArrayCooMap{N} <: Map
  blocksize::NTuple{N,Int}
  blockids::Vector{NTuple{N,Int}}
  ptrs::Array{Int,N}
  function BlockFieldArrayCooMap(blocksize::NTuple{N,Int}, blockids::Vector{NTuple{N,Int}}) where N
    ptrs = fill(-1,blocksize)
    for (p,I) in enumerate(blockids)
      ptrs[I...] = p
      for i in 1:N
        @check blocksize[i] >= I[i]
      end
    end
    new{N}(blocksize,blockids,ptrs)
  end
end

@inline function lazy_map(k::BlockFieldArrayCooMap,T::Type,f::AbstractArray...)
  s = Arrays._common_size(f...)
  N = length(s)
  LazyArray(T,Val(N),Fill(k,s),f...)
end

BlockArrays.blocksize(a::BlockFieldArrayCooMap) = a.blocksize
is_zero_block(a::BlockFieldArrayCooMap,i::Integer) = a.ptrs[i] < 0
is_zero_block(a::BlockFieldArrayCooMap{N},i::Vararg{Integer,N}) where N = a.ptrs[i...] < 0
is_zero_block(a::BlockFieldArrayCooMap{N},i::Vararg{Block,N}) where N = is_zero_block(a,map(Int,i)...)
is_zero_block(a::BlockFieldArrayCooMap,i::Block) = is_zero_block(a,convert(Tuple,i)...)
is_zero_block(a::BlockFieldArrayCooMap,i::CartesianIndex) = is_zero_block(a,Tuple(i)...)

function evaluate!(cache,k::BlockFieldArrayCooMap,axes,block)
  @check map(i->first(blocksize(i)),axes) == k.blocksize "The given axes are not compatible with the given BlockFieldArrayCooMap"
  BlockFieldArrayCoo(axes,k.blockids,block)
end

testitem(a::BlockFieldArrayCoo) = testitem(a.block)

# Specific API

function is_zero_block(a::BlockFieldArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  i != first(a.blockids)
end

# AbstractBlockArray

@inline function BlockArrays.getblock(a::BlockFieldArrayCoo{T,N}, i::Vararg{Integer, N}) where {T,N}
  if i == first(a.blockids)
    a.block
  else
    @notimplemented "Cannot get a zero block from a BlockFieldArrayCoo"
  end
end

# AbstractArray

Base.size(a::BlockFieldArrayCoo) = map(length,Base.axes(a))
Base.axes(a::BlockFieldArrayCoo) = a.axes
Base.IndexStyle(::Type{<:BlockFieldArrayCoo}) = IndexCartesian()

function Base.getindex(a::BlockFieldArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

# Evaluation

function return_cache(a::BlockFieldArrayCoo,x::Point)
  fc = return_cache(a.block,x)
  fx = return_value(a.block,x)
  bsize = map(i->first(blocksize(i)),a.axes)
  k = BlockArrayCooMap(bsize,a.blockids)
  cr = return_cache(k,a.axes,fx)
  (k,fc,cr)
end

@inline function evaluate!(cache,a::BlockFieldArrayCoo,x::Point)
  k,fc,cr = cache
  fx = evaluate!(fc,a.block,x)
  evaluate!(cr,k,a.axes,fx)
end

function return_cache(a::BlockFieldArrayCoo,x::AbstractVector{<:Point})
  fc = return_cache(a.block,x)
  fx = return_value(a.block,x)
  pr = similar_range(first(a.axes),length(x))
  axs = (pr,a.axes...)
  blockids = map(i->(1,i...),a.blockids)
  bsize = map(i->first(blocksize(i)),axs)
  k = BlockArrayCooMap(bsize,blockids)
  cr = return_cache(k,axs,fx)
  (k,fc,axs,cr)
end

@inline function evaluate!(cache,a::BlockFieldArrayCoo,x::AbstractVector{<:Point})
  k,fc,axs,cr = cache
  fx = evaluate!(fc,a.block,x)
  evaluate!(cr,k,axs,fx)
end

# Gradient

function evaluate!(cache,k::Broadcasting{typeof(∇)},a::BlockFieldArrayCoo)
  g = k(a.block)
  BlockFieldArrayCoo(a.axes,a.blockids,g)
end

function evaluate!(cache,k::Broadcasting{typeof(∇∇)},a::BlockFieldArrayCoo)
  g = k(a.block)
  BlockFieldArrayCoo(a.axes,a.blockids,g)
end

# Transpose

function Base.transpose(a::BlockFieldArrayCoo{T,1} where T)
  r = similar_range(first(axes(a)),1)
  axs = (r,axes(a)...)
  blockids = map(i->(1,i...),a.blockids)
  BlockFieldArrayCoo(axs,blockids,transpose(a.block))
end

# Global optimizations

# Evaluation

function lazy_map(
  ::typeof(evaluate),
  a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},
  x::AbstractArray{<:Point})

  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_aix = lazy_map(evaluate,cell_ai,x)
  lazy_map(BlockArrayCooMap(m.blocksize,m.blockids),cell_axs,cell_aix)
end

function lazy_map(
  ::typeof(evaluate),
  a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},
  cell_x::AbstractArray{<:AbstractVector{<:Point}})

  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_aix = lazy_map(evaluate,cell_ai,cell_x)

  function result_axes_on_point_vector(axs,x)
    pr = similar_range(first(axs),length(x))
    axs = (pr,axs...)
  end

  cell_axs_new = lazy_map(result_axes_on_point_vector,cell_axs,cell_x)

  bsize_new = (1,m.blocksize...)
  bids_new = map(i->(1,i...),m.blockids)

  r = lazy_map(BlockArrayCooMap(bsize_new,bids_new),cell_axs_new, cell_aix)
end

# Gradient

function lazy_map(
  k::Broadcasting{typeof(∇)}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})

  m = a.maps.value
  cell_axs, cell_ai = a.args

  cell_gi = lazy_map(k,cell_ai)
  lazy_map(m,cell_axs, cell_gi)
end

function lazy_map(
  k::Broadcasting{typeof(∇∇)}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})

  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_gi = lazy_map(k,cell_ai)
  lazy_map(m,cell_axs,cell_gi)
end

function lazy_map(::typeof(axes),a::LazyArray{<:Fill{Broadcasting{typeof(∇)}}})
  b = a.args[1]
  lazy_map(axes,b)
end

function lazy_map(::typeof(axes),a::LazyArray{<:Fill{Broadcasting{typeof(∇∇)}}})
  b = a.args[1]
  lazy_map(axes,b)
end

# Composition

function lazy_map(
  k::Broadcasting{typeof(∘)}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},b::AbstractArray{<:Field})

  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_gi = lazy_map(k,cell_ai,b)
  lazy_map(m,cell_axs,cell_gi)
end

function lazy_map(::typeof(axes),a::LazyArray{<:Fill{Broadcasting{typeof(∘)}}})
  b = a.args[1]
  lazy_map(axes,b)
end

# Transpose

function lazy_map(
  k::typeof(transpose), a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})

  m = a.maps.value
  cell_axs, cell_ai = a.args

  cell_ait = lazy_map(transpose,cell_ai)

  cell_axs_new = lazy_map(cell_axs) do axs
    @check length(axs) == 1
    r = similar_range(first(axs),1)
    (r,axs...)
  end

  bsize_new = (1,m.blocksize...)
  bids_new = map(i->(1,i...),m.blockids)

  lazy_map(BlockFieldArrayCooMap(bsize_new,bids_new),cell_axs_new, cell_ait)
end

# Operations before evaluating

function lazy_map(
  k::Broadcasting{<:Operation}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})
  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_aif = lazy_map(k,cell_ai)
  lazy_map(m,cell_axs,cell_aif)
end

function lazy_map(
  k::Broadcasting{<:Operation}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},f::AbstractArray{<:Field})
  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_aif = lazy_map(k,cell_ai,f)
  lazy_map(m,cell_axs,cell_aif)
end

function lazy_map(
  k::Broadcasting{<:Operation}, f::AbstractArray{<:Field}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})
  m = a.maps.value
  cell_axs, cell_ai = a.args
  cell_aif = lazy_map(k,f,cell_ai)
  lazy_map(m,cell_axs,cell_aif)
end

# Operations on values

function _get_axes_and_blocks(f)
  f[1], f[2:end]
end
#
# Unary operations
# Assumption: op is a scaling of a
function lazy_map(
  k::BroadcastingFieldOpMap, a::LazyArray{<:Fill{<:BlockArrayCooMap}})
  m = a.maps.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.args)
  cell_blocks_new = map(b->lazy_map(k,b),cell_blocks)
  lazy_map(m,cell_axs,cell_blocks_new...)
end

# Binary test/field or trial/field

# Assumption: op is a scaling of a
function _op_basis_vs_field(k,a,f)
  @check ndims(eltype(a)) > ndims(eltype(f))
  m = a.maps.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.args)
  cell_blocks_new = map(b->lazy_map(k,b,f),cell_blocks)
  lazy_map(m,cell_axs,cell_blocks_new...)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{1}}},
  f::AbstractArray{<:Number})
  _op_basis_vs_field(k,a,f)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
  f::AbstractArray{<:Number})
  _op_basis_vs_field(k,a,f)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
  f::AbstractArray{<:AbstractVector})
  _op_basis_vs_field(k,a,f)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{3}}},
  f::AbstractArray{<:AbstractVector})
  _op_basis_vs_field(k,a,f)
end

# Assumption: op is a scaling of a
function _op_field_vs_basis(k,f,a)
  @check ndims(eltype(a)) > ndims(eltype(f))
  m = a.maps.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.args)
  cell_blocks_new = map(b->lazy_map(k,f,b),cell_blocks)
  lazy_map(m,cell_axs,cell_blocks_new...)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:Number},
  a::LazyArray{<:Fill{BlockArrayCooMap{1}}})
  _op_field_vs_basis(k,f,a)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:Number},
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}})
  _op_field_vs_basis(k,f,a)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:AbstractVector},
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}})
  _op_field_vs_basis(k,f,a)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:AbstractVector},
  a::LazyArray{<:Fill{BlockArrayCooMap{3}}})
  _op_field_vs_basis(k,f,a)
end

# Binary test/test or trial/trial
# Assumption: op is a linear combination of a and b
function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{N}}},
  b::LazyArray{<:Fill{BlockArrayCooMap{N}}}) where N

  ma = a.maps.value
  mb = b.maps.value
  @assert blocksize(ma) == blocksize(mb)

  cell_axs, = a.args

  blocks = []
  blockids = eltype(ma.blockids)[]
  for I in eachblockid(ma)
    if is_nonzero_block(ma,I) || is_nonzero_block(mb,I)
      if is_nonzero_block(ma,I) && is_nonzero_block(mb,I)
        aI = _get_cell_block(ma,a,I)
        bI = _get_cell_block(mb,b,I)
      elseif is_nonzero_block(ma,I)
        aI = _get_cell_block(ma,a,I)
        bI = _lazy_array_of_zero_blocks(_block_type(eltype(b)),cell_axs,I,aI)
      else
        bI = _get_cell_block(mb,b,I)
        aI = _lazy_array_of_zero_blocks(_block_type(eltype(a)),cell_axs,I,bI)
      end
      block = lazy_map(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end

  @assert length(blocks) > 0
  m = BlockArrayCooMap(blocksize(ma),blockids)
  lazy_map(m,cell_axs,blocks...)
end

function _block_type(::Type{<:BlockArrayCoo{T,N,A}}) where {T,N,A}
  A
end

function _get_cell_block(
  m::BlockArrayCooMap,
  a::AbstractArray{<:BlockArrayCoo{T,N,A}},
  b::Block) where {T,N,A}
  cell_axs, cell_blocks = _get_axes_and_blocks(a.args)
  I = convert(Tuple,b)
  p = m.ptrs[I...]
  @assert p > 0
  cell_blocks[p]
end

function _lazy_array_of_zero_blocks(::Type{A},cell_axs,I,aI) where A
  lazy_map(cell_axs) do axs
    laxs = map( local_range, axs, I.n)
    Arrays._zero_block(A,laxs)
  end
end

function _lazy_array_of_zero_blocks(
  ::Type{A},cell_axs,I,aI::LazyArray{<:Fill{BlockArrayCooMap{N}}}) where {A,N}

  cell_axsI, = aI.args
  maI = aI.maps.value
  blocks = []
  for J in eachblockid(maI)
    if is_nonzero_block(maI,J)
      block = _lazy_array_of_zero_blocks(_block_type(A),cell_axsI,J,nothing)
      push!(blocks,block)
    end
  end

  lazy_map(maI,cell_axsI,blocks...)
end

# Binary test/trial
# Assumption: op is a product of a and b

#function lazy_map(
#  k::BroadcastingFieldOpMap,
#  a::LazyArray{<:Fill{BlockArrayCooMap{1}}},
#  b::LazyArray{<:Fill{BlockArrayCooMap{2}}})
#
#end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
  b::LazyArray{<:Fill{BlockArrayCooMap{3}}})

  cell_axs_a, = a.args
  cell_axs_b, = b.args
  cell_axs = lazy_map((a1,a2) -> (a1[1],a1[2],a2[3]),cell_axs_a,cell_axs_b)
  ma = a.maps.value
  mb = b.maps.value

  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(ma.ptrs,2)
  nfield2 = size(mb.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(ma,I1) && is_nonzero_block(mb,I2)
        aI1 = _get_cell_block(ma,a,I1)
        bI2 = _get_cell_block(mb,b,I2)
        block = lazy_map(k,aI1,bI2)
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end

  @assert length(blocks) > 0
  bs = (ma.blocksize[1],ma.blocksize[2],mb.blocksize[3])
  m = BlockArrayCooMap(bs,blockids)
  lazy_map(m,cell_axs,blocks...)
end

#function lazy_map(
#  k::BroadcastingFieldOpMap,
#  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
#  b::LazyArray{<:Fill{BlockArrayCooMap{1}}})
#
#end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{3}}},
  b::LazyArray{<:Fill{BlockArrayCooMap{2}}})

  cell_axs_a, = a.args
  cell_axs_b, = b.args
  cell_axs = lazy_map((a1,a2) -> (a2[1],a2[2],a1[3]),cell_axs_a,cell_axs_b)
  ma = a.maps.value
  mb = b.maps.value

  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(mb.ptrs,2)
  nfield2 = size(ma.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(mb,I1) && is_nonzero_block(ma,I2)
        bI1 = _get_cell_block(mb,b,I1)
        aI2 = _get_cell_block(ma,a,I2)
        block = lazy_map(k,bI1,aI2)
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end

  @assert length(blocks) > 0
  bs = (mb.blocksize[1],mb.blocksize[2],ma.blocksize[3])
  m = BlockArrayCooMap(bs,blockids)
  lazy_map(m,cell_axs,blocks...)
end

# Integration of elem vectors
function lazy_map(k::IntegrationMap,a::LazyArray{<:Fill{BlockArrayCooMap{2}}},w::AbstractArray,j::AbstractArray)
  ma = a.maps.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.args)
  cell_axs_new = lazy_map(a->(a[2],),cell_axs)
  blocks = map(block->lazy_map(k,block,w,j),cell_blocks)
  blockids = [ (ids[2],) for ids in ma.blockids ]
  m = BlockArrayCooMap(ma.blocksize[2:end],blockids)
  lazy_map(m,cell_axs_new,blocks...)
end

## Integration of elem matrices
function lazy_map(k::IntegrationMap,a::LazyArray{<:Fill{BlockArrayCooMap{3}}},w::AbstractArray,j::AbstractArray)
  ma = a.maps.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.args)
  cell_axs_new = lazy_map(a->(a[2],a[3]),cell_axs)
  blocks = map(block->lazy_map(k,block,w,j),cell_blocks)
  blockids = [ (ids[2],ids[3]) for ids in ma.blockids ]
  m = BlockArrayCooMap(ma.blocksize[2:end],blockids)
  lazy_map(m,cell_axs_new,blocks...)
end
#function lazy_map(k::IntMap,f::VectorOfBlockArrayCoo{T,3} where T,w::AbstractArray,j::AbstractArray)
#  ax = lazy_map(a->(a[2],a[3]),f.axes)
#  blocks = map(block->lazy_map(k,block,w,j),f.blocks)
#  blockids = [ (ids[2], ids[3]) for ids in f.blockids ]
#  VectorOfBlockArrayCoo(blocks,blockids,ax)
#end
