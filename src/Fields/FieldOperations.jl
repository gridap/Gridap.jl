
# Operate fields at local level

"""
"""
function operate_fields(op::Function,args...)
  k = FieldOpKernel(op)
  lazy_map_kernel_to_field(k,args...)
end

"""
    field_operation_axes(axs::Tuple...)

Returns the axes after lazy_maping a FieldOpKernel to some arrays with axes `axs...`
"""
function field_operation_axes(axs::Tuple...)
  n = maximum( map(length,axs) )
  rs = []
  for i in 1:n
    r = 1:0
    for ax in axs
      if length(ax) >= i
        if length(r) > 1
          @assert length(ax[i]) == 1 || length(ax[i]) == length(r)
        end
        if length(ax[i]) > length(r)
          r = ax[i]
        end
      end
    end
    push!(rs,r)
  end
  Tuple(rs)
end

function field_operation_metasize(axs::Tuple...)
  n = maximum( map(length,axs) )
  rs = []
  for i in 1:n
    r = 1
    for ax in axs
      if length(ax) >= i
        if ax[i] == (:)
          r = ax[i]
        end
      end
    end
    push!(rs,r)
  end
  Tuple(rs)
end

# Operate fields at global level

"""
"""
function operate_arrays_of_fields(op::Function,args...)
  k = FieldOpKernel(op)
  lazy_map_to_field_array(k,args...)
end

function operate_arrays_of_fields(::Type{T},op::Function,args...) where T
  k = FieldOpKernel(op)
  lazy_map_to_field_array(T,k,args...)
end

# Pre-define some operations
# some of them only make sense for fields, not for arrays of fields.

for op in (:+,:-,:tr, :transpose, :adjoint, :symmetric_part)
  @eval begin

    function ($op)(f::Field)
      operate_fields($op,f)
    end

    #function ($op)(f::AbstractArray{<:Field})
    #  operate_arrays_of_fields($op,f)
    #end

  end
end

for op in (:+,:-,:*,:inner,:outer,:dot)
  @eval begin

    function ($op)(f::Field,g::Field)
      operate_fields($op,f,g)
    end

    #function ($op)(f::AbstractArray{<:Field},g::AbstractArray{<:Field})
    #  operate_arrays_of_fields($op,f,g)
    #end

  end
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

@inline function lazy_map_kernel!(cache,k::FieldOpKernel,args...)
  bk = bcast(k.op)
  lazy_map_kernel!(cache,bk,args...)
end

# Define gradients at local and global level

function lazy_map_kernel_gradient(k::FieldOpKernel,args...)
  @notimplemented "The gradient of the result of operation $(k.op) is not yet implemented."
end

function lazy_map_gradient(k::Valued{<:FieldOpKernel},args...)
  @notimplemented "The gradient of the result of operation $(k.op) is not yet implemented."
end

for op in (:+,:-)
  @eval begin

    # Local level

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},a,b)
      ga = field_gradient(a)
      gb = field_gradient(b)
      lazy_map_kernel_to_field(k,ga,gb)
    end

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},a)
      ga = field_gradient(a)
      lazy_map_kernel_to_field(k,ga)
    end

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},a::Number,b)
      gb = field_gradient(b)
      lazy_map_kernel_to_field(k,gb)
    end

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},b,a::Number)
      gb = field_gradient(b)
      gb
    end

    # Global level

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},a,b)
      ga = field_array_gradient(a)
      gb = field_array_gradient(b)
      lazy_map(k,ga,gb)
    end

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},a)
      ga = field_array_gradient(a)
      lazy_map(k,ga)
    end

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},a::Number,b)
      gb = field_array_gradient(b)
      lazy_map(k,gb)
    end

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},b,a::Number)
      gb = field_array_gradient(b)
      gb
    end

  end
end

for op in (:*,⋅)
  @eval begin

    # Local level

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},a::Number,b)
      gb = field_gradient(b)
      lazy_map_kernel_to_field(k,a,gb)
    end

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},b,a::Number)
      gb = field_gradient(b)
      lazy_map_kernel_to_field(k,gb,a)
    end

    function lazy_map_kernel_gradient(k::FieldOpKernel{typeof($op)},a,b)
      ga = field_gradient(a)
      gb = field_gradient(b)
      f1 = lazy_map_kernel_to_field(k,ga,b)
      f2 = lazy_map_kernel_to_field(k,a,gb)
      operate_fields(+,f1,f2)
    end

    # Global level

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},a::AbstractArray{<:Number},b)
      gb = field_array_gradient(b)
      lazy_map(k,a,gb)
    end

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},b,a::AbstractArray{<:Number})
      gb = field_array_gradient(b)
      lazy_map(k,gb,a)
    end

    function lazy_map_gradient(k::Valued{FieldOpKernel{typeof($op)}},a,b)
      ga = field_array_gradient(a)
      gb = field_array_gradient(b)
      f1 = lazy_map(k,ga,b)
      f2 = lazy_map(k,a,gb)
      operate_arrays_of_fields(+,f1,f2)
    end

  end
end

# Local level

function lazy_map_kernel_gradient(k::FieldOpKernel{typeof(/)},b,a::Number)
  gb = field_gradient(b)
  lazy_map_kernel_to_field(k,gb,a)
end

# Global level

function lazy_map_gradient(k::Valued{FieldOpKernel{typeof(/)}},b,a::Number)
  gb = field_array_gradient(b)
  lazy_map(k,gb,a)
end

# Move the value of a test basis into "trial" state

"""
"""
function trialize_basis(f)
  lazy_map_kernel_to_field(trialize_basis_value,f)
end

"""
"""
function trialize_array_of_bases(af)
  lazy_map_to_field_array(trialize_basis_value,af)
end

@inline function trialize_basis_value(a::AbstractMatrix)
  TrializedMatrix(a)
end

function lazy_map_kernel_gradient(::typeof(trialize_basis_value),a)
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

function lazy_map(f::typeof(trialize_basis_value),a::VectorOfBlockMatrixCoo)
  blocks = map(b->lazy_map(f,b),a.blocks)
  blockids = broadcast(ij->(ij[1],1,ij[2]),a.blockids)
  axs = lazy_map( ax -> (ax[1],blockedrange([1]),ax[2]), a.axes)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

# Unary operations
# Assumption: op is linear wrt a
function lazy_map(k::FieldOpKernel,a::VectorOfBlockArrayCoo)
  blocks = map(b->lazy_map(k,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

# Binary test/field or trial/field
# Assumption: op is linear wrt a
function lazy_map(k::FieldOpKernel,a::VectorOfBlockArrayCoo,f::AbstractArray{<:AbstractVector})
  blocks = map(b->lazy_map(k,b,f), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function lazy_map(k::FieldOpKernel,a::VectorOfBlockArrayCoo,f::AbstractArray{<:Number})
  blocks = map(b->lazy_map(k,b,f), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

# Binary field/test or field/trial
# Assumption: op is linear wrt a
function lazy_map(k::FieldOpKernel,f::AbstractArray{<:AbstractVector},a::VectorOfBlockArrayCoo)
  blocks = map(b->lazy_map(k,f,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

function lazy_map(k::FieldOpKernel,f::AbstractArray{<:Number},a::VectorOfBlockArrayCoo)
  blocks = map(b->lazy_map(k,f,b), a.blocks)
  VectorOfBlockArrayCoo(blocks,a.blockids,a.axes,a.ptrs)
end

# Binary test/test or trial/trial
# Assumption: op is a linear combination of a and b
function lazy_map(
  k::FieldOpKernel,a::VectorOfBlockArrayCoo{Ta,N} where Ta,b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N
  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) || is_nonzero_block(b,I)
      block = lazy_map(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

# Binary + test/test or trial/trial
function lazy_map(
  k::FieldOpKernel{typeof(+)},
  a::VectorOfBlockArrayCoo{Ta,N} where Ta,
  b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N

  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) && is_nonzero_block(b,I)
      block = lazy_map(k,aI,bI)
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
function lazy_map(
  k::FieldOpKernel{typeof(-)},
  a::VectorOfBlockArrayCoo{Ta,N} where Ta,
  b::VectorOfBlockArrayCoo{Tb,N} where Tb) where N

  @assert size(a.ptrs) == size(b.ptrs)
  blocks = []
  blockids = NTuple{N,Int}[]
  for (I,aI) in enumerateblocks(a)
    bI = b[I]
    if is_nonzero_block(a,I) && is_nonzero_block(b,I)
      block = lazy_map(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    elseif is_nonzero_block(a,I)
      block = aI
      push!(blocks,block)
      push!(blockids,I.n)
    elseif is_nonzero_block(b,I)
      block = lazy_map(k,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,a.axes)
end

# Binary test/trial
# Assumption: op is a product of a and b
function lazy_map(
  k::FieldOpKernel,a::VectorOfBlockMatrixCoo,b::VectorOfBlockArrayCoo{Tb,3} where Tb)
  axs = lazy_map( (a1,a2) -> (a1[1],a1[2],a2[3]) ,a.axes,b.axes)
  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(a.ptrs,2)
  nfield2 = size(b.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(a,I1) && is_nonzero_block(b,I2)
        block = lazy_map(k,a[I1],b[I2])
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs)
end

# Binary trial/test
# Assumption: op is a product of a and b
function lazy_map(
  k::FieldOpKernel,a::VectorOfBlockArrayCoo{Tb,3} where Tb,b::VectorOfBlockMatrixCoo)
  axs = lazy_map( (a1,a2) -> (a1[1],a2[2],a1[3]) ,a.axes,b.axes)
  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(b.ptrs,2)
  nfield2 = size(a.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(b,I1) && is_nonzero_block(a,I2)
        block = lazy_map(k,b[I1],a[I2])
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end
  VectorOfBlockArrayCoo(Tuple(blocks),blockids,axs)
end

# General operation
# Assumption only one block per direction. This assumption is OK since it is the only
# efficient situation.
function lazy_map(
  k::FieldOpKernel,
  b::VectorOfBlockArrayCoo,
  c::VectorOfBlockArrayCoo,
  d::Union{VectorOfBlockArrayCoo,AbstractArray{<:AbstractVector}}...)
  a = (b,c,d...)
  I = _general_op_block_id(a...)
  aI = map(_get_first_block,a)
  blocks = (lazy_map(k,aI...),)
  blockids = [I,]
  axs = _general_op_axes(a...)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function _general_op_block_id(a...)
  I = zeros(Int,3)
  I[1] = 1
  for ai in a
    i = _get_first_block_id(ai)
    if I[length(i)] == 0
        I[length(i)] = i[end]
      else
      @notimplementedif I[length(i)] != i[end]
    end
  end
  if I[3] == 0
    (I[1],I[2])
  else
    if I[2] == 0
      I[2] = 1
    end
    Tuple(I)
  end
end

_get_first_block_id(a) = (1,)

function _get_first_block_id(a::VectorOfBlockArrayCoo)
  @notimplementedif length(a.blocks) != 1
  first(a.blockids)
end

_get_first_block(a) = a
_get_first_block(a::VectorOfBlockArrayCoo) = first(a.blocks)

function  _general_op_axes(a...)
  test = nothing
  for ai in a
    if isa(ai,VectorOfBlockArrayCoo{T,2} where T)
      test = ai
      break
    end
  end
  trial = nothing
  for ai in a
    if isa(ai,VectorOfBlockArrayCoo{T,3} where T)
      trial = ai
      break
    end
  end
  if test != nothing && trial != nothing
    axs = lazy_map( (a1,a2) -> (a1[1],a1[2],a2[3]),test.axes,trial.axes)
  elseif test == nothing && trial != nothing
    axs = trial.axes
  elseif test != nothing && trial == nothing
    axs = test.axes
  else
    @unreachable
  end
  axs
end

# Integration of elem vectors
function lazy_map(k::IntKernel,f::VectorOfBlockArrayCoo{T,2} where T,w::AbstractArray,j::AbstractArray)
  ax = lazy_map(a->(a[2],),f.axes)
  blocks = map(block->lazy_map(k,block,w,j),f.blocks)
  blockids = [ (ids[2],) for ids in f.blockids ]
  VectorOfBlockArrayCoo(blocks,blockids,ax)
end

# Integration of elem matrices
function lazy_map(k::IntKernel,f::VectorOfBlockArrayCoo{T,3} where T,w::AbstractArray,j::AbstractArray)
  ax = lazy_map(a->(a[2],a[3]),f.axes)
  blocks = map(block->lazy_map(k,block,w,j),f.blocks)
  blockids = [ (ids[2], ids[3]) for ids in f.blockids ]
  VectorOfBlockArrayCoo(blocks,blockids,ax)
end
