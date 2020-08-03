module BlockArraysCooTests

using BlockArrays
using Gridap
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: IntKernel
using Gridap.TensorValues
import Gridap.Arrays: kernel_cache
import Gridap.Arrays: apply_kernel!
import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Arrays: apply

Gridap.Arrays.array_caches() = ()
Gridap.Arrays.getitems!(::Tuple{},::Tuple{},i) = ()

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

function is_zero_block(a,b::Block)
  i = convert(Tuple,b)
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

Base.size(a::BlockArrayCoo) = map(length,Base.axes(a))

Base.axes(a::BlockArrayCoo) = a.axes

function Base.getindex(a::BlockArrayCoo,i::Integer...)
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

struct FieldOp{T} <: Kernel
  op::T
end

function kernel_cache(k::FieldOp,args...)
  bk = bcast(k.op)
  kernel_cache(bk,args...)
end

@inline function apply_kernel!(cache,k::FieldOp,args...)
  bk = bcast(k.op)
  apply_kernel!(cache,bk,args...)
end

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
    @assert all(map(length,blocks) .==  length(first(blocks)) )

    B = typeof(blocks)
    X = typeof(axes)
    blocks_i = collect(map(testitem,blocks))
    axes_i = testitem(axes)
    t = BlockArrayCoo(blocks_i,blockids,axes_i,ptrs)
    T = typeof(t)
    Z = typeof(zero_blocks)
    new{T,N,B,X,Z}(blocks,blockids,axes,ptrs,zero_blocks)
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

function _compute_zero_blocks_array(blocks,ptrs,axes)
  A = eltype(first(blocks))
  cis = CartesianIndices(ptrs)
  zero_blocks = []
  for ci in cis
    p = ptrs[ci]
    if p<0
      block = apply(axes) do a
        i = Tuple(ci)
        s = map((k,l)->length(k[Block(l)]),a,i)
        b = A(undef,s)
        fill!(b,zero(eltype(A)))
        b
      end
      push!(zero_blocks,block)
    end
  end
  Tuple(zero_blocks)
end

Base.IndexStyle(::Type{<:VectorOfBlockArrayCoo}) = IndexLinear()

Base.size(a::VectorOfBlockArrayCoo) = (length(first(a.blocks)),)

function array_cache(a::VectorOfBlockArrayCoo)
  blocks_i = collect(map(testitem,a.blocks))
  zero_blocks_i = collect(eltype(blocks_i),map(testitem,a.zero_blocks))
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

function is_zero_block(a::VectorOfBlockArrayCoo,i::Integer...)
  p = a.ptrs[i...]
  @assert p != 0
  p < 0
end

function apply(k::FieldOp,a::VectorOfBlockArrayCoo)
  blocks = map(b->apply(k,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function apply(k::FieldOp,a::VectorOfBlockArrayCoo,b::AbstractArray)
  blocks = map(block->apply(k,block,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function apply(k::FieldOp,b::AbstractArray,a::VectorOfBlockArrayCoo)
  blocks = map(block->apply(k,b,block), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function apply(k::FieldOp,a::VectorOfBlockArrayCoo{Ta,N} where Ta, b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N
  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  cis = CartesianIndices(a.ptrs)
  for ci in cis
    I = Block(Tuple(ci))
    if (!is_zero_block(a,I)) || (!is_zero_block(b,I))
      block = apply(k,a[I],b[I])
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

function apply(k::FieldOp{typeof(+)},a::VectorOfBlockArrayCoo{Ta,N} where Ta, b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N
  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  cis = CartesianIndices(a.ptrs)
  for ci in cis
    I = Block(Tuple(ci))
    if (!is_zero_block(a,I)) && (!is_zero_block(b,I))
      block = apply(k,a[I],b[I])
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    elseif !is_zero_block(a,I)
      block = a[I]
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    elseif !is_zero_block(b,I)
      block = b[I]
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

function apply(k::FieldOp{typeof(-)},a::VectorOfBlockArrayCoo{Ta,N} where Ta, b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N
  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  cis = CartesianIndices(a.ptrs)
  for ci in cis
    I = Block(Tuple(ci))
    if (!is_zero_block(a,I)) && (!is_zero_block(b,I))
      block = apply(k,a[I],b[I])
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    elseif !is_zero_block(a,I)
      block = a[I]
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    elseif !is_zero_block(b,I)
      block = apply(k,b[I])
      push!(blocks,block)
      push!(blockids,Tuple(ci))
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

function apply(k::FieldOp,a::VectorOfBlockArrayCoo{Ta,2} where Ta, b::VectorOfBlockArrayCoo{Tb,3} where Tb)
  axs = apply( (a1,a2) -> (a1[1],a1[2],a2[3]) ,a.axes,b.axes)
  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(a.ptrs,2)
  nfield2 = size(b.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if !is_zero_block(a,I1) && !is_zero_block(b,I2)
        block = apply(k,a[I1],b[I2])
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs)
end

function apply(k::FieldOp,a::VectorOfBlockArrayCoo{Ta,3} where Ta, b::VectorOfBlockArrayCoo{Tb,2} where Tb)
  axs = apply( (a1,a2) -> (a1[1],a2[2],a1[3]) ,a.axes,b.axes)
  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(b.ptrs,2)
  nfield2 = size(a.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if !is_zero_block(b,I1) && !is_zero_block(a,I2)
        block = apply(k,b[I1],a[I2])
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs)
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

function apply(k::IntKernel,f::VectorOfBlockArrayCoo{T,2} where T,w::AbstractArray,j::AbstractArray)
  ax = apply(a->(a[2],),f.axes)
  blocks = map(block->apply(k,block,w,j),f.blocks)
  blockids = [ (ids[2],) for ids in f.blockids ]
  VectorOfBlockArrayCoo(blocks,blockids,ax)
end

function apply(k::IntKernel,f::VectorOfBlockArrayCoo{T,3} where T,w::AbstractArray,j::AbstractArray)
  ax = apply(a->(a[2],a[3]),f.axes)
  blocks = map(block->apply(k,block,w,j),f.blocks)
  blockids = [ (ids[2], ids[3]) for ids in f.blockids ]
  VectorOfBlockArrayCoo(blocks,blockids,ax)
end

using Test
using BenchmarkTools
using FillArrays

blocks = [ [1 2; 3 4], [5 6; 7 8; 9 10] ]
blockids = [(1,1),(2,1)]
ax = (blockedrange([2,3]), blockedrange([2,4]))
ptrs = _compute_ptrs(blockids,ax)
zero_blocks = _compute_zero_blocks(eltype(blocks),ptrs,ax)

a = BlockArrayCoo(blocks,blockids,ax)
@test a[Block(1),Block(1)] === blocks[1]
@test a[Block(1,1)] === blocks[1]
@test a[Block(1,2)] === a.zero_blocks[1]

b21 = zeros(3,2)
getblock!(b21,a,Block(2,1))
@test b21 == a[Block(2,1)]

b12 = ones(2,4)
getblock!(b12,a,Block(1,2))
@test b12 == a[Block(1,2)]

blocks = [ [1,2,3] ]
blockids = [(2,)]
ax = (blockedrange([2,3]),)
a = BlockArrayCoo(blocks,blockids,ax)

l = 10
b11 = [  i*[1 2; 3 4]  for i in 1:l ]
b21 = [  i*[5 6; 7 8; 9 10] for i in 1:l]
blocks = (b11,b21)
blockids = [(1,1),(2,1)]
ax = Fill((blockedrange([2,3]), blockedrange([2,4])),l)
ptrs = _compute_ptrs(blockids,first(ax))
zero_blocks = _compute_zero_blocks_array(blocks,ptrs,ax)

al = VectorOfBlockArrayCoo(blocks,blockids,ax)
@test al[Block(1,1)] === b11

rl = [ BlockArrayCoo([b11i,b21i],blockids,axi,ptrs)  for (b11i,b21i,axi) in zip(b11,b21,ax) ]
test_array(al,rl)

npoin = 3
ndofs = 4
nfiel = 2
ncell = 10

# single-field
cf = [rand(npoin) for cell in 1:ncell]
cbv = [rand(npoin,ndofs) for cell in 1:ncell]
cbu = trialize_array_of_matrices(cbv)

cj = [ fill(TensorValue(1,0,0,2),npoin) for cell in 1:ncell ]

# random operations between cell fields
f(f1,f2) = f1*f2+f1
ca = apply(FieldOp(f),cf,cf)
cr = [ f.(cfi,cfi)  for cfi in cf]
test_array(ca,cr)

# Linear operations on cell basis (test)
f(u,u2,f) = u +f*u2
ca = apply(FieldOp(f),cbv,cbv,cf)
cr = [ f.(cbvi,cbvi,cfi)  for (cbvi,cfi) in zip(cbv,cf)]
@test size(ca[1]) == (npoin,ndofs)
test_array(ca,cr)

# Linear operations on cell basis (trial)
f(u,u2,f) = u +f*u2
ca = apply(FieldOp(f),cbu,cbu,cf)
cr = [ f.(cbui,cbui,cfi)  for (cbui,cfi) in zip(cbu,cf)]
@test size(ca[1]) == (npoin,1,ndofs)
test_array(ca,cr)

# Bi-linear operations on test/trial cell bases
f(v,v2,u,u2,f) = (v+v2)*(u+f*u2)
ca = apply(FieldOp(f),cbv,cbv,cbu,cbu,cf)
cr = [ f.(cbvi,cbvi,cbui,cbui,cfi)  for (cbvi,cbui,cfi) in zip(cbv,cbu,cf)]
@test size(ca[1]) == (npoin,ndofs,ndofs)
test_array(ca,cr)

# multi-field
cf1 = [rand(npoin) for cell in 1:ncell]
cf2 = [rand(npoin) for cell in 1:ncell]
b1 = [rand(npoin,ndofs) for cell in 1:ncell]
b2 = [rand(npoin,ndofs) for cell in 1:ncell]
axs = Fill( (blockedrange([npoin]),blockedrange([ndofs,ndofs])),ncell)
cbv1 = VectorOfBlockArrayCoo((b1,),[(1,1)],axs)
cbv2 = VectorOfBlockArrayCoo((b2,),[(1,2)],axs)
cbu1 = trialize_array_of_matrices(cbv1)
cbu2 = trialize_array_of_matrices(cbv2)

# Unary operations on cell basis (test)
f(u) = VectorValue(u,u)
ca = apply(FieldOp(f),cbv1)

f(a,b) = VectorValue(a,b)
ca = apply(FieldOp(f),cbv1,cf1)
ca = apply(FieldOp(f),cf1,cbv1)

# Unary operations on cell basis (trial)
f(u) = VectorValue(u,u)
ca = apply(FieldOp(f),cbu1)

f(a,b) = VectorValue(a,b)
ca = apply(FieldOp(f),cbu1,cf1)
ca = apply(FieldOp(f),cf1,cbu1)

# Linear operations on cell basis (test)
f(u,u2) = u + 2*u2
ca = apply(FieldOp(f),cbv1,cbv2)

ca = apply(FieldOp(+),cbv1,cbv2)
@test ca[Block(1,1)] === cbv1[Block(1,1)]
@test ca[Block(1,2)] === cbv2[Block(1,2)]

ca = apply(FieldOp(-),cbv1,cbv2)
@test ca[Block(1,1)] === cbv1[Block(1,1)]

# Linear operations on cell basis (trial)
f(u,u2) = u + 2*u2
ca = apply(FieldOp(f),cbu1,cbu2)

ca = apply(FieldOp(+),cbu1,cbu2)
@test ca[Block(1,1,1)] === cbu1[Block(1,1,1)]
@test ca[Block(1,1,2)] === cbu2[Block(1,1,2)]

ca = apply(FieldOp(-),cbu1,cbu2)
@test ca[Block(1,1,1)] === cbu1[Block(1,1,1)]

# Bi-linear operations on cell basis (test/trial)
f(v,u) = 2*v*u
ca = apply(FieldOp(f),cbv1,cbu2)
@test ca.blockids == [(1,1,2)]

# Bi-linear operations on cell basis (trial/test)
ca = apply(FieldOp(f),cbu1,cbv2)
@test ca.blockids == [(1,2,1)]

# Integration of vectors
ca = apply(IntKernel(),cbv1,cf1,cj)
display(ca.axes)
display(ca[1])

ca = apply(IntKernel(),apply(FieldOp(*),cbv1,cbu2),cf1,cj)
display(ca[1])

end # module
