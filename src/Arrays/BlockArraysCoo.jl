
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

function BlockArrayCoo(
  blocks::Vector{A},ax::NTuple{N},ptrs::Array{Int,N},zero_blocks::Vector{A}) where {A,N}
  cis = CartesianIndices(ptrs)
  blockids = NTuple{N,Int}[]
  for ci in cis
    p = ptrs[ci]
    if p>0
      push!(blockids,Tuple(ci))
    end
  end
  BlockArrayCoo(blocks,blockids,ax,ptrs,zero_blocks)
end

function BlockArrayCoo(allblocks::AbstractArray{A,N},ax::NTuple{N},mask::AbstractArray{Bool,N}=fill(true,size(allblocks))) where {A,N}
  blocks = A[]
  zero_blocks = A[]
  ptrs = zeros(Int,size(allblocks))
  pp = 1
  pn = 1
  cis = CartesianIndices(allblocks)
  for ci in cis
    if mask[ci]
      push!(blocks,allblocks[ci])
      ptrs[ci] = pp
      pp += 1
    else
      push!(zero_blocks,allblocks[ci])
      ptrs[ci] = -pn
      pn += 1
    end
  end
  BlockArrayCoo(blocks,ax,ptrs,zero_blocks)
end

#  axs = _compute_axes_from_allblocks(allblocks)
#
#function _compute_axes_from_allblocks(allblocks)
#  N = ndims(allblocks)
#  s = size(allblocks)
#  axtmp = [ zeros(Int,s[n]) for n in 1:N]
#  cis = CartesianIndices(allblocks)
#  for ci in cis
#    block = allblocks[ci]
#    sb = size(block)
#    for (n,cin) in enumerate(Tuple(ci))
#      axtmp[n][cin] = sb[n]
#    end
#  end
#  axs = map(blockedrange,Tuple(axtmp))
#  axs
#end

const BlockMatrixCoo = BlockArrayCoo{T,2} where T
const BlockVectorCoo = BlockArrayCoo{T,1} where T

zeros_like(a::AbstractArray) = zeros_like(a,axes(a))

function zeros_like(a::BlockArrayCoo,axs::Tuple)
  blocks = eltype(a.blocks)[]
  blockids = eltype(a.blockids)[]
  BlockArrayCoo(blocks,blockids,axs)
end

function zeros_like(a::AbstractArray,axs::Tuple)
  T = eltype(a)
  b = similar(a,T,axs)
  fill!(b,zero(T))
  b
end

function zeros_like(::Type{<:BlockArrayCoo{T,N,A}},axs::Tuple) where {T,N,A}
  blocks = A[]
  blockids = NTuple{N,Int}[]
  BlockArrayCoo(blocks,blockids,axs)
end

function zeros_like(::Type{A},axs::Tuple) where A <: AbstractArray
  a = similar(A,axs)
  fill!(a,zero(eltype(a)))
  a
end

Base.similar(a::BlockArrayCoo) = similar(a,eltype(a),axes(a))

Base.similar(a::BlockArrayCoo,axs::Tuple) = similar(a,eltype(a),axs)

function Base.similar(a::BlockArrayCoo{S,N} where S,::Type{T}, axs::NTuple{N,<:BlockedUnitRange}) where {T,N}
  _similar_block_array_coo(a,T,axs)
end

function _similar_block_array_coo(a,::Type{T},axs) where T
  @notimplementedif length(a.blocks) == 0
  A = typeof(similar(first(a.blocks),T))
  blocks = A[]
  for p in 1:length(a.blocks)
    I = a.blockids[p]
    laxs = map( local_range, axs, I)
    block = similar(a.blocks[p],T,laxs)
    push!(blocks,block)
  end
  BlockArrayCoo(blocks,a.blockids,axs,a.ptrs)
end

function local_range(a::BlockedUnitRange,k::Integer)
  n = length(a[Block(k)])
  Base.OneTo(n)
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
      block = zeros_like(A,laxs)
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
  blocks = eachblockindex(a)
  zip(blocks,eachblock(a))
end

function eachblockindex(a)
  cis = CartesianIndices(blocksize(a))
  map(ci->Block(Tuple(ci)),cis)
end

Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))

Base.axes(a::BlockArrayCoo) = a.axes

function Base.getindex(a::BlockArrayCoo{T,N} where T,i::Vararg{Integer,N}) where N
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

function Base.getindex(a::BlockVectorCoo,i::Integer,j::Integer...)
  @assert all( j .== 1)
  s = map(findblockindex,a.axes,(i,))
  ai = a[s...]
  ai
end

function Base.setindex!(a::BlockArrayCoo{T,N} where T,v,i::Vararg{Integer,N}) where N
  s = map(findblockindex,a.axes,i)
  s = map(findblockindex,a.axes,(i,))
  I = Block(map(i->i.I[1],s)...)
  α = CartesianIndex(map(BlockArrays.blockindex,s))
  a[I][α] = v
end

function Base.setindex!(a::BlockVectorCoo,v,i::Integer,j::Integer...)
  @assert all( j .== 1)
  s = map(findblockindex,a.axes,(i,))
  I = Block(map(i->i.I[1],s)...)
  α = CartesianIndex(map(BlockArrays.blockindex,s))
  a[I][α] = v
end

function Base.:+(a::BlockArrayCoo{Ta,N},b::BlockArrayCoo{Tb,N}) where {Ta,Tb,N}
  @assert axes(a) == axes(b)
  @assert blockaxes(a) == blockaxes(b)
  I1 = first(eachblockindex(a))
  A = typeof(a[I1]+b[I1])
  blocks = A[]
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) || is_nonzero_block(b,I)
      block = aI + bI
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  BlockArrayCoo(blocks,blockids,a.axes)
end

function Base.:-(a::BlockArrayCoo{Ta,N},b::BlockArrayCoo{Tb,N}) where {Ta,Tb,N}
  @assert axes(a) == axes(b)
  @assert blockaxes(a) == blockaxes(b)
  I1 = first(eachblockindex(a))
  A = typeof(a[I1]-b[I1])
  blocks = A[]
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) || is_nonzero_block(b,I)
      block = aI - bI
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  BlockArrayCoo(blocks,blockids,a.axes)
end

function Base.:*(a::BlockArrayCoo,b::Number)
  *(b,a)
end

function Base.:*(a::Number,b::BlockArrayCoo)
  f(block)= a*block
  blocks = f.(b.blocks)
  BlockArrayCoo(blocks,b.blockids,b.axes,b.ptrs,b.zero_blocks)
end

function Base.:*(a::BlockMatrixCoo,b::BlockVectorCoo)
  @assert blocksize(a,2) == blocksize(b,1)
  A = typeof(a[Block(1,1)]*b[Block(1)])
  blocks = A[]
  blockids = Tuple{Int}[]
  for i in 1:blocksize(a,1)
    block = zeros_like(a[Block(i,1)]*b[Block(1)])
    for j in 1:blocksize(a,2)
      if is_nonzero_block(a,Block(i,j)) && is_nonzero_block(b,Block(j))
        block += a[Block(i,j)]*b[Block(j)]
      end
    end
    push!(blocks,block)
    push!(blockids,(i,))
  end
  axs = (axes(a)[1],)
  BlockArrayCoo(blocks,blockids,axs)
end

function Base.:*(a::BlockMatrixCoo,b::BlockMatrixCoo)
  @assert blocksize(a,2) == blocksize(b,1)
  A = typeof(a[Block(1,1)]*b[Block(1,1)])
  blocks = A[]
  blockids = Tuple{Int,Int}[]
  for i in 1:blocksize(a,1)
    for j in 1:blocksize(b,2)
      block = zeros_like(a[Block(i,1)]*b[Block(1,j)])
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
  BlockArrayCoo(blocks,blockids,axs)
end

function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo)
  fill!(c,zero(eltype(c)))
  mul!(c,a,b,1,0)
end

function LinearAlgebra.mul!(c::BlockMatrixCoo,a::BlockMatrixCoo,b::BlockMatrixCoo)
  fill!(c,zero(eltype(c)))
  mul!(c,a,b,1,0)
end

function scale_entries!(c::BlockArrayCoo,β)
  for block in c.blocks
    scale_entries!(block,β)
  end
  c
end

function LinearAlgebra.mul!(c::BlockVectorCoo,a::BlockMatrixCoo,b::BlockVectorCoo,α::Number,β::Number)
  for I in 1:blocksize(a,1)
    cI = c[Block(I)]
    if is_nonzero_block(c,Block(I))
      scale_entries!(cI,β)
    end
    for J in 1:blocksize(a,2)
      if is_nonzero_block(a,Block(I,J)) && is_nonzero_block(b,Block(J))
        @assert is_nonzero_block(c,Block(I))
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
          @assert is_nonzero_block(c,Block(I,J))
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

function Base.fill!(a::BlockArrayCoo,v)
  @notimplementedif v != zero(v)
  for b in a.blocks
    fill!(b,v)
  end
  a
end

function LinearAlgebra.transpose(a::BlockMatrixCoo)
  blocks = transpose.(a.blocks)
  zero_blocks = transpose.(a.zero_blocks)
  blockids = [ (j,i) for (i,j) in a.blockids ]
  ax,ay = axes(a)
  axs = (ay,ax)
  ptrs = collect(transpose(a.ptrs))
  BlockArrayCoo(blocks,blockids,axs,ptrs,zero_blocks)
end

function Base.copy(a::BlockArrayCoo)
  blocks = copy.(a.blocks)
  blockids = copy(a.blockids)
  axes = copy.(a.axes)
  ptrs = copy(a.ptrs)
  zero_blocks = copy.(a.zero_blocks)
  BlockArrayCoo(blocks,blockids,axes,ptrs,zero_blocks)
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
    for I in eachblockindex(a)
      if is_nonzero_block(a,I)
        copyto!(a[I],b[I])
      end
    end
  end
 a
end

# For blocks of blocks

struct TwoLevelBlockedUnitRange{CS} <: AbstractUnitRange{Int}
  global_range::BlockedUnitRange{CS}
  local_ranges::Vector{BlockedUnitRange{CS}}
end

function TwoLevelBlockedUnitRange(local_ranges::Vector{<:BlockedUnitRange})
  global_range = blockedrange(length.(local_ranges))
  TwoLevelBlockedUnitRange(global_range,local_ranges)
end

function BlockArrays.blockedrange(local_ranges::Vector{<:BlockedUnitRange})
  TwoLevelBlockedUnitRange(local_ranges)
end

function Base.similar(a::BlockArrayCoo{S,N} where S,::Type{T}, axs::NTuple{N,<:TwoLevelBlockedUnitRange}) where {T,N}
  _similar_block_array_coo(a,T,axs)
end

function Base.similar(::Type{T},axs::NTuple{N,<:TwoLevelBlockedUnitRange}) where {N,T<:Array}
  similar(T,map(length,axs))
end

local_range(a::TwoLevelBlockedUnitRange,k::Integer) = a.local_ranges[k]

Base.first(a::TwoLevelBlockedUnitRange) = first(a.global_range)
Base.last(a::TwoLevelBlockedUnitRange) = last(a.global_range)
BlockArrays.blockaxes(a::TwoLevelBlockedUnitRange) = blockaxes(a.global_range)
BlockArrays.blocklasts(a::TwoLevelBlockedUnitRange) = blocklasts(a.global_range)
BlockArrays.findblock(a::TwoLevelBlockedUnitRange,k::Integer) = findblock(a.global_range,k)
Base.getindex(a::TwoLevelBlockedUnitRange,i::Integer) = a.global_range[i]
Base.getindex(a::TwoLevelBlockedUnitRange,i::Block{1}) = a.global_range[i]
Base.axes(a::TwoLevelBlockedUnitRange) = axes(a.global_range)
Base.Broadcast.axistype(a::T, b::T) where T<:TwoLevelBlockedUnitRange =  Base.Broadcast.axistype(a.global_range,b.global_range)
Base.Broadcast.axistype(a::TwoLevelBlockedUnitRange, b::TwoLevelBlockedUnitRange) = Base.Broadcast.axistype(a.global_range,b.global_range)
Base.Broadcast.axistype(a::TwoLevelBlockedUnitRange, b) = Base.Broadcast.axistype(a.global_range,b)
Base.Broadcast.axistype(a, b::TwoLevelBlockedUnitRange) =Base.Broadcast.axistype(a,b.global_range)
Base.print_matrix_row(
  io::IO,
  X::TwoLevelBlockedUnitRange,
  A::Vector,
  i::Integer,
  cols::AbstractVector,
  sep::AbstractString) = print_matrix_row(io,X.global_range,A,i,cols,sep)
Base.show(io::IO, mimetype::MIME"text/plain", a::TwoLevelBlockedUnitRange) = show(io,mimetype,a.global_range)

