module BlockArraysCooTests

using BlockArrays
using Gridap
using Gridap.Arrays

struct BlockArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  blocks::Vector{A}
  blockids::Vector{NTuple{N,Int}}
  axes::X
  ptrs::Array{Int,N}
  function BlockArrayCoo(
    blocks::Vector{<:AbstractArray{T,N}},
    blockids::Vector{NTuple{N,Int}},
    axes::NTuple{N},
    ptrs::Array{Int,N}=_compute_ptrs(blockids,axes)) where {T,N}

    A = eltype(blocks)
    X = typeof(axes)
    new{T,N,A,X}(blocks,blockids,axes,ptrs)
  end
end

function _compute_ptrs(blockids,axes)
  s = map(i->first(blocksize(i)),axes)
  ptrs = zeros(Int,s)
  for (i,c) in enumerate(blockids)
    ptrs[c...] = i
  end
  ptrs
end

function BlockArrays.getblock(a::BlockArrayCoo{T,N,A} where {T,N},i::Integer...) where A
  p = a.ptrs[i...]
  if p>0
    a.blocks[p]
  else
    s = map((k,l)->length(k[Block(l)]),a.axes,i)
    block = A(undef,s)
    fill!(block,zero(eltype(A)))
    block
  end
end

function BlockArrays.getblock!(block,a::BlockArrayCoo{T,N,A} where {T,N},i::Integer...) where A
  p = a.ptrs[i...]
  if p>0
    block .= a.blocks[p]
  else
    fill!(block,zero(eltype(A)))
    block
  end
end

Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))

Base.axes(a::BlockArrayCoo) = a.axes

function Base.getindex(a::BlockArrayCoo,i::Integer...)
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

struct VectorOfBlockArrayCoo{T,B,N,X} <: AbstractVector{T}
  blocks::B
  blockids::Vector{NTuple{N,Int}}
  axes::X
  ptrs::Array{Int,N}
  function VectorOfBlockArrayCoo(
    blocks::Tuple,
    blockids::Vector{NTuple{N,Int}},
    axes::AbstractArray{<:NTuple{N}},
    ptrs::Array{Int,N} = _compute_ptrs(blockids,testitem(axes))) where N

    msg = "Trying to build a VectorOfBlockArrayCoo with repeated blocks"
    @assert _no_has_repeaded_blocks(blockids,ptrs) msg
    @assert all(map(length,blocks) .==  length(first(blocks)) )

    B = typeof(blocks)
    X = typeof(axes)
    blocks_i = collect(map(testitem,blocks))
    axes_i = testitem(axes)
    t = BlockArrayCoo(blocks_i,blockids,axes_i,ptrs)
    T = typeof(t)
    new{T,B,N,X}(blocks,blockids,axes,ptrs)
  end
end

function _no_has_repeaded_blocks(blockids::Vector{NTuple{N,Int}},ptrs) where N
  maxblocks = size(ptrs)
  touched = zeros(Int,maxblocks)
  for c in blockids
    touched[c...] += 1
  end
  all( touched .<= 1 )
end

Base.IndexStyle(::Type{<:VectorOfBlockArrayCoo}) = IndexLinear()

Base.size(a::VectorOfBlockArrayCoo) = (length(first(a.blocks)),)

function Gridap.Arrays.array_cache(a::VectorOfBlockArrayCoo)
  blocks_i = collect(map(testitem,a.blocks))
  ca = array_cache(a.axes)
  cb = array_caches(a.blocks...)
  (blocks_i,ca,cb)
end

@inline function Gridap.Arrays.getindex!(cache,a::VectorOfBlockArrayCoo,i::Integer)
  blocks_i, ca, cb = cache
  axes_i = getindex!(ca,a.axes,i)
  blocks_i .= getitems!(cb,a.blocks,i)
  BlockArrayCoo(blocks_i,a.blockids,axes_i,a.ptrs)
end

function Base.getindex(a::VectorOfBlockArrayCoo,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end


struct TrializedMatrix{T,A} <: AbstractArray{T,3}
  matrix::A
  @inline function TrializedMatrix(matrix::AbstractMatrix{T}) where T
    A = typeof(matrix)
    new{T,A}(matrix)
  end
end

TrializedMatrix{T,A}(u::UndefInitializer,s::Tuple) where {T,A} = TrializedMatrix(A(u,(s[1],s[3])))

Base.size(a::TrializedMatrix) = (size(a.matrix,1),1,size(a.matrix,2))

Base.IndexStyle(::Type{<:TrializedMatrix{T,A}}) where {T,A} = IndexStyle(A)

@inline Base.getindex(a::TrializedMatrix,i::Integer,j::Integer,k::Integer) = a.matrix[i,k]

@inline Base.getindex(a::TrializedMatrix,i::Integer) = a.matrix[i]

@inline Base.setindex!(a::TrializedMatrix,v,i::Integer,j::Integer,k::Integer) = (a.matrix[i,k] = v)

@inline Base.setindex!(a::TrializedMatrix,v,i::Integer) = (a.matrix[i] = v)

@inline function trialize_matrix(a::AbstractMatrix)
  TrializedMatrix(a)
end

function trialize_array_of_matrices(a)
  apply(trialize_matrix,a)
end

function trialize_array_of_matrices(a::VectorOfBlockArrayCoo)
  blocks = map(trialize_array_of_matrices,a.blocks)
  blockids = broadcast(ij->(ij[1],1,ij[2]),a.blockids)
  axs = apply( ax -> (ax[1],blockedrange([1]),ax[2]), a.axes)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

using Test
using BenchmarkTools
using FillArrays

blocks = [ [1 2; 3 4], [5 6; 7 8; 9 10] ]
blockids = [(1,1),(2,1)]
axes = (blockedrange([2,3]), blockedrange([2,4]))
ptrs = _compute_ptrs(blockids,axes)

a = BlockArrayCoo(blocks,blockids,axes)
@test a[Block(1),Block(1)] === blocks[1]
@test a[Block(1,1)] === blocks[1]

b21 = zeros(3,2)
getblock!(b21,a,Block(2,1))
@test b21 == a[Block(2,1)]

b12 = ones(2,4)
getblock!(b12,a,Block(1,2))
@test b12 == a[Block(1,2)]

l = 10
b11 = [  i*[1 2; 3 4]  for i in 1:l ]
b21 = [  i*[5 6; 7 8; 9 10] for i in 1:l]
blocks = (b11,b21)
blockids = [(1,1),(2,1)]
axes = Fill((blockedrange([2,3]), blockedrange([2,4])),l)
ptrs = _compute_ptrs(blockids,first(axes))

al = VectorOfBlockArrayCoo(blocks,blockids,axes)

rl = [ BlockArrayCoo([b11i,b21i],blockids,axi,ptrs)  for (b11i,b21i,axi) in zip(b11,b21,axes) ]
test_array(al,rl)

npoin = 3
ndofs = 4
nfiel = 2
ncell = 10

# single-field
cf = [rand(npoin) for cell in 1:ncell]
cb = [rand(npoin,ndofs) for cell in 1:ncell]

# multi-field
cf1 = [rand(npoin) for cell in 1:ncell]
cf2 = [rand(npoin) for cell in 1:ncell]
b1 = [rand(npoin,ndofs) for cell in 1:ncell]
b2 = [rand(npoin,ndofs) for cell in 1:ncell]
axs = Fill( (blockedrange([npoin]),blockedrange([ndofs,ndofs])),ncell)
cb1 = VectorOfBlockArrayCoo((b1,),[(1,1)],axs)
cb2 = VectorOfBlockArrayCoo((b2,),[(1,2)],axs)

#cell = 2
#display(cb[cell])
#display(cb1[cell])
#display(cb2[cell])

a = rand(3,4)
b = trialize_matrix(a)
b2 = reshape(a,(3,1,4))
test_array(b,b2)

cbt = trialize_array_of_matrices(cb)

cb1t = trialize_array_of_matrices(cb1)
cache = array_cache(cb1t)
@btime getindex!($cache,$cb1t,1)
#display(cb1t)

cb2t = trialize_array_of_matrices(cb2)
#display(cb2t)





end # module
