# Blocked axes

function append_ranges(ranges::AbstractVector{<:AbstractUnitRange})
  blockedrange(map(length,ranges))
end

function append_ranges(ranges::AbstractVector{<:BlockedUnitRange})
  MultiLevelBlockedUnitRange(ranges)
end

"""
Implementation of a multilevel blocked range
It is a generalization of BlockedUnitRange that stores more information
about the local ranges.
"""
struct MultiLevelBlockedUnitRange{L} <: AbstractUnitRange{Int}
  local_ranges::Vector{L}
  global_range::BlockedUnitRange{Vector{Int}}
  function MultiLevelBlockedUnitRange(local_ranges::Vector{L}) where L<:AbstractUnitRange
    local_lengths = map(length,local_ranges)
    global_range = blockedrange(local_lengths)
    new{L}(local_ranges,global_range)
  end
end

Base.first(a::MultiLevelBlockedUnitRange) = first(a.global_range)
Base.last(a::MultiLevelBlockedUnitRange) = last(a.global_range)
BlockArrays.blockaxes(a::MultiLevelBlockedUnitRange) = blockaxes(a.global_range)
BlockArrays.blocklasts(a::MultiLevelBlockedUnitRange) = blocklasts(a.global_range)
BlockArrays.findblock(a::MultiLevelBlockedUnitRange,k::Integer) = findblock(a.global_range,k)
Base.getindex(a::MultiLevelBlockedUnitRange,i::Integer) = a.global_range[i]
Base.getindex(a::MultiLevelBlockedUnitRange,i::Block{1}) = a.global_range[i]
Base.getindex(a::MultiLevelBlockedUnitRange, i::BlockRange{1}) = a.global_range[i]
Base.axes(a::MultiLevelBlockedUnitRange) = axes(a.global_range)
Base.Broadcast.axistype(a::T, b::T) where T<:MultiLevelBlockedUnitRange =  Base.Broadcast.axistype(a.global_range,b.global_range)
Base.Broadcast.axistype(a::MultiLevelBlockedUnitRange, b::MultiLevelBlockedUnitRange) = Base.Broadcast.axistype(a.global_range,b.global_range)
Base.Broadcast.axistype(a::MultiLevelBlockedUnitRange, b) = Base.Broadcast.axistype(a.global_range,b)
Base.Broadcast.axistype(a, b::MultiLevelBlockedUnitRange) =Base.Broadcast.axistype(a,b.global_range)
#Base.print_matrix_row(
#  io::IO,
#  X::MultiLevelBlockedUnitRange,
#  A::Vector,
#  i::Integer,
#  cols::AbstractVector,
#  sep::AbstractString) = print_matrix_row(io,X.global_range,A,i,cols,sep)

local_range(a::AbstractUnitRange,k::Integer) = Base.OneTo(length(a[Block(k)]))
local_range(a::MultiLevelBlockedUnitRange,k::Integer) = a.local_ranges[k]

function similar_range(r::Base.OneTo,n::Integer)
  Base.OneTo(Int(n))
end

function similar_range(r::BlockedUnitRange,n::Integer)
  blockedrange([n])
end

function similar_range(r::MultiLevelBlockedUnitRange,n::Integer)
  r = similar_range(first(r.local_ranges),n)
  append_ranges([r])
end

"""
Check if the full multi-level block structure is the same.
This is in contrast to BlockArrays.blockisequal that only checks
the top level block structure.
"""
@inline blocks_equal(a::AbstractUnitRange,b::AbstractUnitRange) = BlockArrays.blockisequal(a,b)
blocks_equal(a::AbstractUnitRange,b::MultiLevelBlockedUnitRange) = false
blocks_equal(a::MultiLevelBlockedUnitRange,b::AbstractUnitRange) = false
function blocks_equal(a::MultiLevelBlockedUnitRange,b::MultiLevelBlockedUnitRange)
  if a === b
    return true
  end
  r = blocks_equal(a.global_range,b.global_range)
  if r == false
    return false
  end
  la = length(a.local_ranges)
  lb = length(b.local_ranges)
  if la!=lb
    return false
  else
    for i in 1:la
      @inbounds ra = a.local_ranges[i]
      @inbounds rb = b.local_ranges[i]
      r = r && blocks_equal(ra,rb)
    end
    return r
  end
end


"""
Check if the number ob blocks in all levels is the same.
"""
@inline function num_blocks_equal(a::AbstractUnitRange,b::AbstractUnitRange)
  sa = BlockArrays.blocksize(a)
  sb = BlockArrays.blocksize(b)
  sa == sb
end
num_blocks_equal(a::AbstractUnitRange,b::MultiLevelBlockedUnitRange) = false
num_blocks_equal(a::MultiLevelBlockedUnitRange,b::AbstractUnitRange) = false
function num_blocks_equal(a::MultiLevelBlockedUnitRange,b::MultiLevelBlockedUnitRange)
  if a === b
    return true
  end
  if length(a.local_ranges) != length(b.local_ranges)
    return false
  end
  for i in 1:length(a.local_ranges)
    if ! num_blocks_equal(a.local_ranges[i],b.local_ranges[i])
      return false
    end
  end
  true
end

blocks_equal(a::Tuple,b::Tuple) = all(map(blocks_equal,a,b))
num_blocks_equal(a::Tuple,b::Tuple) = all(map(num_blocks_equal,a,b))

"""
Block array in which some of the blocks are known to be zero
"""
struct BlockArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  axes::X
  blockids::Vector{NTuple{N,Int}}
  blocks::Vector{A}
  ptrs::Array{Int,N}
  zero_blocks::Vector{A}
  function BlockArrayCoo(
    axes::NTuple{N},
    blockids::Vector{NTuple{N,Int}},
    blocks::Vector{A},
    ptrs::Array{Int,N},
    zero_blocks::Vector{A}) where {T,N,A<:AbstractArray{T,N}}

    @assert length(blockids) == length(blocks)
    @check length(unique(blockids)) == length(blockids) "We cannot built a BlockArrayCoo from repeated blocks"
    @check _valid_block_sizes(axes,ptrs,blocks) "The given blocks do not match with the given axes"

    X = typeof(axes)
    new{T,N,A,X}(axes,blockids,blocks,ptrs,zero_blocks)
  end
end

function _valid_block_sizes(_axes,ptrs,blocks)
  s = map(i->first(blocksize(i)),_axes)
  cis = CartesianIndices(s)
  for ci in cis
    p = ptrs[ci]
    if p>0
      I = Tuple(ci)
      laxs = map( local_range, _axes, I)
      if !blocks_equal(axes(blocks[p]),laxs)
        return false
      end
    end
  end
  true
end

# Minimal constructor

function BlockArrayCoo(
  axes::NTuple{N},
  blockids::Vector{NTuple{N,Int}},
  blocks::Vector{A}) where {T,N,A<:AbstractArray{T,N}}

  ptrs = _compute_ptrs(blockids,axes)
  zero_blocks = _compute_zero_blocks(A,ptrs,axes)
  BlockArrayCoo(axes,blockids,blocks,ptrs,zero_blocks)
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

function _compute_zero_blocks(::Type{A},ptrs,axs) where A
  cis = CartesianIndices(ptrs)
  zero_blocks = A[]
  for ci in cis
    p = ptrs[ci]
    if p<0
      I = Tuple(ci)
      laxs = map( local_range, axs, I)
      block = _zero_block(A,laxs)
      push!(zero_blocks,block)
    end
  end
  zero_blocks
end

function _zero_block(::Type{A},axs::Tuple) where A <: AbstractArray
  a = similar(A,axs)
  fill!(a,zero(eltype(a)))
  a
end

function _zero_block(::Type{<:Transpose{T,A}},axs::Tuple) where {T,A}
  a = similar(A,(axs[2],))
  fill!(a,zero(eltype(a)))
  Transpose(a)
end

function _zero_block(::Type{<:BlockArrayCoo{T,N,A}},axs::Tuple) where {T,N,A}
  blocks = A[]
  blockids = NTuple{N,Int}[]
  BlockArrayCoo(axs,blockids,blocks)
end

# Minimal constructor (for lazy_map)

struct BlockArrayCooMap{N} <: Map
  blocksize::NTuple{N,Int}
  blockids::Vector{NTuple{N,Int}}
  ptrs::Array{Int,N}
  function BlockArrayCooMap(blocksize::NTuple{N,Int}, blockids::Vector{NTuple{N,Int}}) where N
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

@inline function lazy_map(k::BlockArrayCooMap,T::Type,f::AbstractArray...)
  s = _common_size(f...)
  N = length(s)
  LazyArray(T,Val(N),Fill(k,s),f...)
end

#function BlockArrayCoo(
#  axes::NTuple{N},
#  blockids::Vector{NTuple{N,Int}},
#  blocks::A...) where {T,N,A<:AbstractArray{T,N}}
#
#  BlockArrayCoo(axes,blockids,collect(blocks))
#end

function return_cache(
  k::BlockArrayCooMap,
  axes::NTuple{N},
  blocks::A...) where {T,N,A<:AbstractArray{T,N}}

  @check map(i->first(blocksize(i)),axes) == k.blocksize "The given axes are not compatible with the given BlockArrayCooMap"

  if _valid_block_sizes(axes,_compute_ptrs(k.blockids,axes),blocks)
    goodblocks = collect(blocks)
  else
    goodblocks = A[]
    for (i,I) in enumerate(k.blockids)
      laxs = map( local_range, axes, I)
      block = similar(typeof(blocks[i]),laxs)
      push!(goodblocks,block)
    end
  end

  r = BlockArrayCoo(axes,k.blockids,goodblocks)
  CachedArray(r)
end

@inline function evaluate!(
  cache,
  ::BlockArrayCooMap,
  axes::NTuple{N},
  blocks::A...) where {T,N,A<:AbstractArray{T,N}}
                           
  setaxes!(cache,axes)
  r = cache.array
  copyto!(r.blocks,blocks)
  r
end

BlockArrays.blocksize(a::BlockArrayCooMap) = a.blocksize
is_zero_block(a::BlockArrayCooMap,i::Integer) = a.ptrs[i] < 0
is_zero_block(a::BlockArrayCooMap{N},i::Vararg{Integer,N}) where N = a.ptrs[i...] < 0
is_zero_block(a::BlockArrayCooMap{N},i::Vararg{Block,N}) where N = is_zero_block(a,map(Int,i)...)
is_zero_block(a::BlockArrayCooMap,i::Block) = is_zero_block(a,convert(Tuple,i)...)
is_zero_block(a::BlockArrayCooMap,i::CartesianIndex) = is_zero_block(a,Tuple(i)...)

# Specific API

@inline is_nonzero_block(args...) = ! is_zero_block(args...)

function is_zero_block(a::AbstractArray,b::Block)
  i = convert(Tuple,b)
  is_zero_block(a,i...)
end

function is_zero_block(a::AbstractArray{T,N},b::Vararg{Block,N}) where {T,N}
  i = map(Int,b)
  is_zero_block(a,i...)
end

function is_zero_block(a::AbstractArray,b::CartesianIndex)
  i = Tuple(b)
  is_zero_block(a,i...)
end

function is_zero_block(a::BlockArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  #@boundscheck BlockArrays.blockcheckbounds(a, i...)
  p = a.ptrs[i...]
  @check p != 0
  p < 0
end

function is_zero_block(a::AbstractArray{T,N},i::Vararg{Integer,N}) where {T,N}
  @boundscheck BlockArrays.blockcheckbounds(a, i...)
  false
end

function eachblockid(a)
  cis = CartesianIndices(blocksize(a))
  map(ci->Block(Tuple(ci)),cis)
end

BlockArrays.eachblock(a::BlockArrayCoo) = ( a[I]  for I in eachblockid(a) )
enumerateblocks(a) = zip(eachblockid(a),eachblock(a))

# AbstractBlockArray

@inline function BlockArrays.getblock(a::BlockArrayCoo{T,N}, block::Vararg{Integer, N}) where {T,N}
  #@boundscheck BlockArrays.blockcheckbounds(a, block...)
  p = a.ptrs[block...]
  if p>0
    a.blocks[p]
  else
    a.zero_blocks[-p]
  end
end

function BlockArrays.getblock!(c,a::BlockArrayCoo{T,N}, block::Vararg{Integer, N}) where {T,N}
  b = a[Block(block...)]
  copy!(c,b)
  c
end

@inline function BlockArrays.setblock!(a::BlockArrayCoo{T, N}, v, block::Vararg{Integer, N}) where {T,N}
  #@boundscheck BlockArrays.blockcheckbounds(a, block...)
  p = a.ptrs[block...]
  if p>0
    a.blocks[p] = v
  else
    @unreachable "It is not possible to set a zero block in BlockArrayCoo"
  end
  a
end

# AbstractArray

Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))
Base.axes(a::BlockArrayCoo) = a.axes
Base.IndexStyle(::Type{<:BlockArrayCoo}) = IndexCartesian()

function Base.getindex(a::BlockArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

function Base.setindex!(a::BlockArrayCoo{T,N},v,i::Vararg{Integer,N}) where {T,N}
  s = map(findblockindex,a.axes,i)
  I = Block(map(i->i.I[1],s)...)
  α = CartesianIndex(map(BlockArrays.blockindex,s))
  a[I][α] = v
end

# Fast copy

function Base.copy(a::BlockArrayCoo)
  blocks = copy.(a.blocks)
  blockids = copy(a.blockids)
  axes = copy.(a.axes)
  ptrs = copy(a.ptrs)
  zero_blocks = copy.(a.zero_blocks)
  BlockArrayCoo(axes,blockids,blocks,ptrs,zero_blocks)
end

function Base.copy!(a::BlockArrayCoo,b::BlockArrayCoo)
  copyto!(a,b)
end

function Base.copyto!(a::BlockArrayCoo,b::BlockArrayCoo)
  if a.ptrs === b.ptrs || a.ptrs == b.ptrs
    for p in 1:length(a.blocks)
      copyto!(a.blocks[p],b.blocks[p])
    end
  else
    for I in eachblockid(a)
      if is_nonzero_block(a,I)
        copyto!(a[I],b[I])
      end
    end
  end
 a
end

# Similar
# We customize similar only for matching axis type
# Otherwise we do not have enough info to create a BlockArrayCoo

function Base.similar(a::BlockArrayCoo)
  similar(a,eltype(a),axes(a))
end

function Base.similar(a::BlockArrayCoo,::Type{T}) where T
  similar(a,T,axes(a))
end

function Base.similar(a::BlockArrayCoo{S,N,A,X},::Type{T}, _axes::X) where {S,A,T,N,X}
  _similar_block_array_coo(a,T,_axes)
end

function Base.similar(
  a::BlockArrayCoo{S,N,A,X},::Type{T}, _axes::X) where {S,A,T,N,X<:Tuple{Vararg{BlockedUnitRange}}}
  _similar_block_array_coo(a,T,_axes)
end

function _similar_block_array_coo(a,::Type{T},_axes) where T
  if num_blocks_equal(axes(a),_axes)
    _similar_block_array_coo_preserving(a,T,_axes)
  else
    _similar_block_array_coo_non_preserving(a,T,_axes)
  end
end

# similar preserving zero block structure
function _similar_block_array_coo_preserving(a,::Type{T},_axes) where T
  ai = first(eachblock(a))
  A = typeof(similar(ai,T))
  blocks = A[]
  for i in 1:length(a.blocks)
    ai = a.blocks[i]
    I = a.blockids[i]
    laxs = map( local_range, _axes, I)
    block = similar(ai,T,laxs)
    push!(blocks,block)
  end
  blockids = copy(a.blockids)
  BlockArrayCoo(_axes,blockids,blocks)
end

# similar without zero block structure (all blocks are assumed to be non-zero)
function _similar_block_array_coo_non_preserving(a::BlockArrayCoo{S,N,A},::Type{T},_axes) where {S,N,A,T}
  s = map(i->first(blocksize(i)),_axes)
  cis = CartesianIndices(s)
  blocks = A[]
  blockids = NTuple{N,Int}[]
  for ci in cis
    I = Tuple(ci)
    laxs = map( local_range, _axes, I)
    block = similar(similar(A,laxs),T,laxs)
    push!(blocks,block)
    push!(blockids,I)
  end
  BlockArrayCoo(_axes,blockids,blocks)
end

# In this case we can not preserve the zero block structure neither (all blocks are assumed to be non-zero).
function Base.similar(::Type{<:BlockArrayCoo{S,N,A,X}},_axes::X) where {S,N,A,X}
  s = map(i->first(blocksize(i)),_axes)
  cis = CartesianIndices(s)
  blocks = A[]
  blockids = NTuple{N,Int}[]
  for ci in cis
    I = Tuple(ci)
    laxs = map( local_range, _axes, I)
    block = similar(A,laxs)
    push!(blocks,block)
    push!(blockids,I)
  end
  BlockArrayCoo(_axes,blockids,blocks)
end

# More API

function Base.fill!(a::BlockArrayCoo,v)
  @notimplementedif v != zero(v)
  for b in a.blocks
    fill!(b,v)
  end
  a
end

function fill_entries!(a::BlockArrayCoo,v)
  for b in a.blocks
    fill_entries!(b,v)
  end
  a
end

function scale_entries!(c::BlockArrayCoo,β)
  for block in c.blocks
    scale_entries!(block,β)
  end
  c
end

Base.zero(a::BlockArrayCoo) = _zero_block(typeof(a),axes(a))

# + and -

for op in (:+,:-)
  @eval begin
    function Base.$op(a::BlockArrayCoo{Ta,N},b::BlockArrayCoo{Tb,N}) where {Ta,Tb,N}
      @check blocks_equal(axes(a),axes(b))
      I1 = first(eachblockid(a))
      A = typeof($op(a[I1],b[I1]))
      blocks = A[]
      blockids = NTuple{N,Int}[]
      for (I,aI) in enumerateblocks(a)
        bI = b[I]
        if is_nonzero_block(a,I) || is_nonzero_block(b,I)
          block = $op(aI,bI)
          push!(blocks,block)
          push!(blockids,I.n)
        end
      end
      BlockArrayCoo(a.axes,blockids,blocks)
    end
  end
end

Base.:+(a::BlockArrayCoo) = a

function Base.:-(a::BlockArrayCoo)
  b = copy(a)
  scale_entries!(b,-1)
  b
end

# Products

function Base.:*(a::BlockArrayCoo,b::Number)
  *(b,a)
end

function Base.:*(a::Number,b::BlockArrayCoo)
  f(block)= a*block
  blocks = f.(b.blocks)
  BlockArrayCoo(b.axes,b.blockids,blocks,b.ptrs,f.(b.zero_blocks))
end

const BlockMatrixCoo = BlockArrayCoo{T,2} where T
const BlockVectorCoo = BlockArrayCoo{T,1} where T

function Base.:*(a::BlockMatrixCoo,b::BlockVectorCoo)
  @assert blocksize(a,2) == blocksize(b,1)
  A = typeof(a[Block(1,1)]*b[Block(1)])
  blocks = A[]
  blockids = Tuple{Int}[]
  for i in 1:blocksize(a,1)
    block = zero(a[Block(i,1)]*b[Block(1)])
    for j in 1:blocksize(a,2)
      if is_nonzero_block(a,Block(i,j)) && is_nonzero_block(b,Block(j))
        block += a[Block(i,j)]*b[Block(j)]
      end
    end
    push!(blocks,block)
    push!(blockids,(i,))
  end
  axs = (axes(a)[1],)
  BlockArrayCoo(axs,blockids,blocks)
end

function Base.:*(a::BlockMatrixCoo,b::BlockMatrixCoo)
  @assert blocksize(a,2) == blocksize(b,1)
  A = typeof(a[Block(1,1)]*b[Block(1,1)])
  blocks = A[]
  blockids = Tuple{Int,Int}[]
  for i in 1:blocksize(a,1)
    for j in 1:blocksize(b,2)
      block = zero(a[Block(i,1)]*b[Block(1,j)])
      for k in 1:blocksize(a,2)
        if is_nonzero_block(a,Block(i,k)) && is_nonzero_block(b,Block(k,j))
          block += a[Block(i,k)]*b[Block(k,j)]
        end
      end
      push!(blocks,block)
      push!(blockids,(i,j))
    end
  end
  axs = (axes(a)[1],axes(b)[2])
  BlockArrayCoo(axs,blockids,blocks)
end

function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo)
  fill_entries!(c,zero(eltype(c)))
  mul!(c,a,b,1,0)
end

function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo)
  fill_entries!(c,zero(eltype(c)))
  mul!(c,a,b,1,0)
end

function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo,α::Number,β::Number)
  for I in 1:blocksize(a,1)
    cI = c[Block(I)]
    if is_nonzero_block(c,Block(I))
      scale_entries!(cI,β)
    end
    for J in 1:blocksize(a,2)
      if is_nonzero_block(a,Block(I,J)) && is_nonzero_block(b,Block(J))
        @check is_nonzero_block(c,Block(I))
        aIJ = a[Block(I,J)]
        bJ = b[Block(J)]
        mul!(cI,aIJ,bJ,α,1)
      end
    end
  end
  c
end

function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo,α::Number,β::Number)
  for I in 1:blocksize(a,1)
    for J in 1:blocksize(b,2)
      cIJ = c[Block(I,J)]
      if is_nonzero_block(c,Block(I,J))
        scale_entries!(cIJ,β)
      end
      for K in 1:blocksize(a,2)
        if is_nonzero_block(a,Block(I,K)) && is_nonzero_block(b,Block(K,J))
          @check is_nonzero_block(c,Block(I,J))
          cIJ = c[Block(I,J)]
          aIK = a[Block(I,K)]
          bKJ = b[Block(K,J)]
          mymul!(cIJ,aIK,bKJ,α) # Hack to avoid allocations
          #mul!(cIJ,aIK,bKJ,α,1) # Why this leads to memory allocations?? Specially with Ints
        end
      end
    end
  end
  c
end

# Hack to avoid allocations
@inline function mymul!(cIJ,aIK,bKJ,α)
  @boundscheck begin
    @assert size(cIJ,1) == size(aIK,1)
    @assert size(cIJ,2) == size(bKJ,2)
    @assert size(aIK,2) == size(bKJ,1)
  end
  for i in 1:size(cIJ,1)
    for j in 1:size(bKJ,2)
      for k in 1:size(bKJ,1)
        @inbounds cIJ[i,j] += α*aIK[i,k]*bKJ[k,j]
      end
    end
  end
end

# Hack to avoid allocations 
@inline function mymul!(cIJ::BlockMatrixCoo,aIK::BlockMatrixCoo,bKJ::BlockMatrixCoo,α)
  mul!(cIJ,aIK,bKJ,α,1)
end

function LinearAlgebra.transpose(a::BlockMatrixCoo)
  blocks = transpose.(a.blocks)
  zero_blocks = transpose.(a.zero_blocks)
  blockids = [ (j,i) for (i,j) in a.blockids ]
  ax,ay = axes(a)
  axs = (ay,ax)
  ptrs = collect(transpose(a.ptrs))
  BlockArrayCoo(axs,blockids,blocks,ptrs,zero_blocks)
end




























#  function BlockArrayCoo(
#    blocks::Vector{A},
#    blockids::Vector{NTuple{N,Int}},
#    axes::NTuple{N},
#    ptrs::Array{Int,N}=_compute_ptrs(blockids,axes),
#    zero_blocks::Vector{A}=_compute_zero_blocks(A,ptrs,axes)) where {T,N,A<:AbstractArray{T,N}}
#
#    X = typeof(axes)
#    new{T,N,A,X}(blocks,blockids,axes,ptrs,zero_blocks)
#  end
#
#function BlockArrayCoo(
#  blocks::Vector{A},ax::NTuple{N},ptrs::Array{Int,N},zero_blocks::Vector{A}) where {A,N}
#  cis = CartesianIndices(ptrs)
#  blockids = NTuple{N,Int}[]
#  for ci in cis
#    p = ptrs[ci]
#    if p>0
#      push!(blockids,Tuple(ci))
#    end
#  end
#  BlockArrayCoo(blocks,blockids,ax,ptrs,zero_blocks)
#end
#
#function BlockArrayCoo(allblocks::AbstractArray{A,N},ax::NTuple{N},mask::AbstractArray{Bool,N}=fill(true,size(allblocks))) where {A,N}
#  blocks = A[]
#  zero_blocks = A[]
#  ptrs = zeros(Int,size(allblocks))
#  pp = 1
#  pn = 1
#  cis = CartesianIndices(allblocks)
#  for ci in cis
#    if mask[ci]
#      push!(blocks,allblocks[ci])
#      ptrs[ci] = pp
#      pp += 1
#    else
#      push!(zero_blocks,allblocks[ci])
#      ptrs[ci] = -pn
#      pn += 1
#    end
#  end
#  BlockArrayCoo(blocks,ax,ptrs,zero_blocks)
#end
#
##function BlockArrayCoo()
##end
#
##  axs = _compute_axes_from_allblocks(allblocks)
##
##function _compute_axes_from_allblocks(allblocks)
##  N = ndims(allblocks)
##  s = size(allblocks)
##  axtmp = [ zeros(Int,s[n]) for n in 1:N]
##  cis = CartesianIndices(allblocks)
##  for ci in cis
##    block = allblocks[ci]
##    sb = size(block)
##    for (n,cin) in enumerate(Tuple(ci))
##      axtmp[n][cin] = sb[n]
##    end
##  end
##  axs = map(blockedrange,Tuple(axtmp))
##  axs
##end
#
#
#zeros_like(a::AbstractArray) = zeros_like(a,axes(a))
#
#function zeros_like(a::BlockArrayCoo,axs::Tuple)
#  blocks = eltype(a.blocks)[]
#  blockids = eltype(a.blockids)[]
#  BlockArrayCoo(blocks,blockids,axs)
#end
#
#function zeros_like(a::AbstractArray,axs::Tuple)
#  T = eltype(a)
#  b = similar(a,T,axs)
#  fill!(b,zero(T))
#  b
#end
#
#function zeros_like(::Type{<:BlockArrayCoo{T,N,A}},axs::Tuple) where {T,N,A}
#  blocks = A[]
#  blockids = NTuple{N,Int}[]
#  BlockArrayCoo(blocks,blockids,axs)
#end
#
#function zeros_like(::Type{A},axs::Tuple) where A <: AbstractArray
#  a = similar(A,axs)
#  fill!(a,zero(eltype(a)))
#  a
#end
#
#
#function local_range(a::BlockedUnitRange,k::Integer)
#  n = length(a[Block(k)])
#  Base.OneTo(n)
#end
#
#
#
#
#function BlockArrays.getblock!(block,a::BlockArrayCoo,i::Integer...)
#  p = a.ptrs[i...]
#  if p>0
#    block .= a.blocks[p]
#  else
#    block .= a.zero_blocks[-p]
#  end
#end
#
#function BlockArrays.eachblock(a::BlockArrayCoo)
#  cis = CartesianIndices(blocksize(a))
#  blocks = map(ci->Block(Tuple(ci)),cis)
#  ( a[block] for block in blocks  )
#end
#
#
#Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))
#
#Base.axes(a::BlockArrayCoo) = a.axes
#
#function Base.getindex(a::BlockArrayCoo{T,N} where T,i::Vararg{Integer,N}) where N
#  s = map(findblockindex,a.axes,i)
#  ai = a[s...]
#  ai
#end
#
#function Base.getindex(a::BlockVectorCoo,i::Integer,j::Integer...)
#  @assert all( j .== 1)
#  s = map(findblockindex,a.axes,(i,))
#  ai = a[s...]
#  ai
#end
#
#function Base.setindex!(a::BlockArrayCoo{T,N} where T,v,i::Vararg{Integer,N}) where N
#  s = map(findblockindex,a.axes,i)
#  s = map(findblockindex,a.axes,(i,))
#  I = Block(map(i->i.I[1],s)...)
#  α = CartesianIndex(map(BlockArrays.blockindex,s))
#  a[I][α] = v
#end
#
#function Base.setindex!(a::BlockVectorCoo,v,i::Integer,j::Integer...)
#  @assert all( j .== 1)
#  s = map(findblockindex,a.axes,(i,))
#  I = Block(map(i->i.I[1],s)...)
#  α = CartesianIndex(map(BlockArrays.blockindex,s))
#  a[I][α] = v
#end
#
#function Base.fill!(a::BlockArrayCoo,v)
#  @notimplementedif v != zero(v)
#  for b in a.blocks
#    fill!(b,v)
#  end
#  a
#end
#
#
#
## For blocks of blocks
#
#
#struct TwoLevelBlockedUnitRange{CS} <: AbstractUnitRange{Int}
#  global_range::BlockedUnitRange{CS}
#  local_ranges::Vector{BlockedUnitRange{CS}}
#end
#
#function TwoLevelBlockedUnitRange(local_ranges::Vector{<:BlockedUnitRange})
#  global_range = blockedrange(length.(local_ranges))
#  TwoLevelBlockedUnitRange(global_range,local_ranges)
#end
#
#function BlockArrays.blockedrange(local_ranges::Vector{<:BlockedUnitRange})
#  TwoLevelBlockedUnitRange(local_ranges)
#end
#
#
#local_range(a::TwoLevelBlockedUnitRange,k::Integer) = a.local_ranges[k]
#
#Base.first(a::TwoLevelBlockedUnitRange) = first(a.global_range)
#Base.last(a::TwoLevelBlockedUnitRange) = last(a.global_range)
#BlockArrays.blockaxes(a::TwoLevelBlockedUnitRange) = blockaxes(a.global_range)
#BlockArrays.blocklasts(a::TwoLevelBlockedUnitRange) = blocklasts(a.global_range)
#BlockArrays.findblock(a::TwoLevelBlockedUnitRange,k::Integer) = findblock(a.global_range,k)
#Base.getindex(a::TwoLevelBlockedUnitRange,i::Integer) = a.global_range[i]
#Base.getindex(a::TwoLevelBlockedUnitRange,i::Block{1}) = a.global_range[i]
#Base.axes(a::TwoLevelBlockedUnitRange) = axes(a.global_range)
#Base.Broadcast.axistype(a::T, b::T) where T<:TwoLevelBlockedUnitRange =  Base.Broadcast.axistype(a.global_range,b.global_range)
#Base.Broadcast.axistype(a::TwoLevelBlockedUnitRange, b::TwoLevelBlockedUnitRange) = Base.Broadcast.axistype(a.global_range,b.global_range)
#Base.Broadcast.axistype(a::TwoLevelBlockedUnitRange, b) = Base.Broadcast.axistype(a.global_range,b)
#Base.Broadcast.axistype(a, b::TwoLevelBlockedUnitRange) =Base.Broadcast.axistype(a,b.global_range)
#Base.print_matrix_row(
#  io::IO,
#  X::TwoLevelBlockedUnitRange,
#  A::Vector,
#  i::Integer,
#  cols::AbstractVector,
#  sep::AbstractString) = print_matrix_row(io,X.global_range,A,i,cols,sep)
#Base.show(io::IO, mimetype::MIME"text/plain", a::TwoLevelBlockedUnitRange) = show(io,mimetype,a.global_range)
#
