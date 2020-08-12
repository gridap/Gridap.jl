
# Operate fields at local level

"""
"""
function operate_fields(op::Function,args...)
  k = FieldOpKernel(op)
  apply_kernel_to_field(k,args...)
end

# Define gradients at local level

function apply_kernel_gradient(k::FieldOpKernel,args...)
  @notimplemented "The gradient of the result of operation $(k.op) is not yet implemented."
end

for op in (:+,:-)
  @eval begin

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},a,b)
      ga = field_gradient(a)
      gb = field_gradient(b)
      apply_kernel_to_field(k,ga,gb)
    end

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},a)
      ga = field_gradient(a)
      apply_kernel_to_field(k,ga)
    end

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},a::Number,b)
      gb = field_gradient(b)
      apply_kernel_to_field(k,gb)
    end

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},b,a::Number)
      gb = field_gradient(b)
      apply_kernel_to_field(k,gb)
    end

  end
end

for op in (:*,⋅)
  @eval begin

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},a::Number,b)
      gb = field_gradient(b)
      apply_kernel_to_field(k,a,gb)
    end

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},b,a::Number)
      gb = field_gradient(b)
      apply_kernel_to_field(k,gb,a)
    end

    function apply_kernel_gradient(k::FieldOpKernel{typeof($op)},a,b)
      ga = field_gradient(a)
      gb = field_gradient(b)
      f1 = apply_kernel_to_field(k,ga,b)
      f2 = apply_kernel_to_field(k,a,gb)
      operate_fields(+,f1,f2)
    end

  end
end

function apply_kernel_gradient(k::FieldOpKernel{typeof(/)},b,a::Number)
  gb = field_gradient(b)
  apply_kernel_to_field(k,gb,a)
end

# Operate fields at global level

"""
"""
function operate_arrays_of_fields(op::Function,args...)
  k = FieldOpKernel(op)
  apply_to_field_array(k,args...)
end

function operate_arrays_of_fields(::Type{T},op::Function,args...) where T
  k = FieldOpKernel(op)
  apply_to_field_array(T,k,args...)
end

# Operations between field values

struct FieldOpKernel{T} <: Kernel
  op::T
end

# FieldOpKernel does broadcast by default. It will work always assuming that the trial bases have shape (np,1,ndof)
# but perhaps inefficient for blocked matrices.
# In any case, optimizations for block matrices will be done at the global level (for all cells)
# instead of at the cell level in this kernel.
# In other words, we can assume that this kernel receives standard non-blocked arrays in practice.

function kernel_cache(k::FieldOpKernel,args...)
  bk = bcast(k.op)
  kernel_cache(bk,args...)
end

@inline function apply_kernel!(cache,k::FieldOpKernel,args...)
  bk = bcast(k.op)
  apply_kernel!(cache,bk,args...)
end

# Move the value of a test basis into "trial" state

"""
"""
function trialize_basis(f)
  apply_kernel_to_field(trialize_basis_value,f)
end

"""
"""
function trialize_array_of_bases(af)
  apply_to_field_array(trialize_basis_value,af)
end

@inline function trialize_basis_value(a::AbstractMatrix)
  TrializedMatrix(a)
end

function apply_kernel_gradient(::typeof(trialize_basis_value),a)
  g = field_gradient(a)
  trialize_basis(g)
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

# Optimizations for block matrices at global level (for all cells)

function apply(f::typeof(trialize_basis_value),a::VectorOfBlockMatrixCoo)
  blocks = map(b->apply(f,b),a.blocks)
  blockids = broadcast(ij->(ij[1],1,ij[2]),a.blockids)
  axs = apply( ax -> (ax[1],blockedrange([1]),ax[2]), a.axes)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

# Unary operations
# Assumption: op is linear wrt a
function apply(k::FieldOpKernel,a::VectorOfBlockArrayCoo)
  blocks = map(b->apply(k,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

# Binary test/field or trial/field
# Assumption: op is linear wrt a
function apply(k::FieldOpKernel,a::VectorOfBlockArrayCoo,f::AbstractArray{<:AbstractVector})
  blocks = map(b->apply(k,b,f), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function apply(k::FieldOpKernel,a::VectorOfBlockArrayCoo,f::AbstractArray{<:Number})
  blocks = map(b->apply(k,b,f), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

# Binary field/test or field/trial
# Assumption: op is linear wrt a
function apply(k::FieldOpKernel,f::AbstractArray{<:AbstractVector},a::VectorOfBlockArrayCoo)
  blocks = map(b->apply(k,f,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function apply(k::FieldOpKernel,f::AbstractArray{<:Number},a::VectorOfBlockArrayCoo)
  blocks = map(b->apply(k,f,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

# Binary test/test or trial/trial
# Assumption: op is a linear combination of a and b
function apply(
  k::FieldOpKernel,a::VectorOfBlockArrayCoo{Ta,N} where Ta,b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N
  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) || is_nonzero_block(b,I)
      block = apply(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

# Binary + test/test or trial/trial
function apply(
  k::FieldOpKernel{typeof(+)},
  a::VectorOfBlockArrayCoo{Ta,N} where Ta,
  b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N

  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) && is_nonzero_block(b,I)
      block = apply(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    elseif is_nonzero_block(a,I)
      block = aI
      push!(blocks,block)
      push!(blockids,I.n)
    elseif is_nonzero_block(b,I)
      block = bI
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

# Binary - test/test or trial/trial
function apply(
  k::FieldOpKernel{typeof(-)},
  a::VectorOfBlockArrayCoo{Ta,N} where Ta,
  b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N

  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) && is_nonzero_block(b,I)
      block = apply(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    elseif is_nonzero_block(a,I)
      block = aI
      push!(blocks,block)
      push!(blockids,I.n)
    elseif is_nonzero_block(b,I)
      block = apply(k,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

# Binary test/trial
# Assumption: op is a product of a and b
function apply(
  k::FieldOpKernel,a::VectorOfBlockMatrixCoo,b::VectorOfBlockArrayCoo{Tb,3} where Tb)
  axs = apply( (a1,a2) -> (a1[1],a1[2],a2[3]) ,a.axes,b.axes)
  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(a.ptrs,2)
  nfield2 = size(b.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(a,I1) && is_nonzero_block(b,I2)
        block = apply(k,a[I1],b[I2])
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs)
end

# Binary trial/test
# Assumption: op is a product of a and b
function apply(
  k::FieldOpKernel,a::VectorOfBlockArrayCoo{Tb,3} where Tb,b::VectorOfBlockMatrixCoo)
  axs = apply( (a1,a2) -> (a1[1],a2[2],a1[3]) ,a.axes,b.axes)
  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(b.ptrs,2)
  nfield2 = size(a.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(b,I1) && is_nonzero_block(a,I2)
        block = apply(k,b[I1],a[I2])
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs)
end

# General operation
# TODO


# Integration of elem vectors
function apply(k::IntKernel,f::VectorOfBlockArrayCoo{T,2} where T,w::AbstractArray,j::AbstractArray)
  ax = apply(a->(a[2],),f.axes)
  blocks = map(block->apply(k,block,w,j),f.blocks)
  blockids = [ (ids[2],) for ids in f.blockids ]
  VectorOfBlockArrayCoo(blocks,blockids,ax)
end

# Integration of elem matrices
function apply(k::IntKernel,f::VectorOfBlockArrayCoo{T,3} where T,w::AbstractArray,j::AbstractArray)
  ax = apply(a->(a[2],a[3]),f.axes)
  blocks = map(block->apply(k,block,w,j),f.blocks)
  blockids = [ (ids[2], ids[3]) for ids in f.blockids ]
  VectorOfBlockArrayCoo(blocks,blockids,ax)
end

