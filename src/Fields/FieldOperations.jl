"""
    field_operation(op::Function,a)
    field_operation(op::Function,a,b)
"""
function field_operation(op::Function,a,b)
  k = FieldBinOp(op)
  apply_kernel_to_field(k,a,b)
end

function field_operation(op::Function,a)
  k = bcast(op)
  apply_kernel_to_field(k,a)
end

"""
    field_array_operation(op::Function,a)
    field_array_operation(op::Function,a,b)
    field_array_operation(::Type{T},op::Function,a) where T
    field_array_operation(::Type{T},op::Function,a,b) where T
"""
function field_array_operation(op::Function,a,b)
  k = FieldBinOp(op)
  apply_to_field_array(k,a,b)
end

function field_array_operation(op::Function,a)
  k = bcast(op)
  apply_to_field_array(k,a)
end

function field_array_operation(::Type{T},op::Function,a) where T
  k = bcast(op)
  apply_to_field_array(T,k,a)
end

function field_array_operation(::Type{T},op::Function,a,b) where T
  k = FieldBinOp(op)
  apply_to_field_array(T,k,a,b)
end

# Unary operations on fields and arrays of fields

for op in (:+,:-,:tr, :transpose, :adjoint, :symmetric_part)
  @eval begin

    function ($op)(f::Field)
      field_operation($op,f)
    end

    function ($op)(f::AbstractArray{<:Field})
      field_array_operation($op,f)
    end

  end

end

# Binary operations on fields and arrays of fields

for op in (:+,:-,:*,:inner,:outer)
  @eval begin

    function ($op)(f::Field, g::Field)
      field_operation($op,f,g)
    end

    function ($op)(f::AbstractArray{<:Field},g::AbstractArray{<:Field})
      field_array_operation($op,f,g)
    end

  end

end

# Helpers

struct FieldBinOp{F<:Function} <: Kernel
  op::F
end

# Vector vs Vector

function kernel_cache(
  k::FieldBinOp,a::AbstractVector,b::AbstractVector)
  na = length(a)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(k.op,Ta,Tb)
  r = zeros(T,na)
  CachedArray(r)
end

function kernel_testitem!(
  c,k::FieldBinOp,a::AbstractVector,b::AbstractVector)
  if _valid_checks_vecvec(a,b)
    apply_kernel!(c,k,a,b)
  else
    c.array
  end
end

@inline function apply_kernel!(
  c,k::FieldBinOp,a::AbstractVector,b::AbstractVector)
  _field_bin_op_checks_vecvec(a,b)
  na = length(a)
  setsize!(c,(na,))
  r = c.array
  for p in eachindex(a)
    @inbounds r[p] = k.op(a[p],b[p])
  end
  r
end

function _field_bin_op_checks_vecvec(a,b)
  @assert _valid_checks_vecvec(a,b) "Binary operation between fields: vector vs vector size mismatch."
end

function _valid_checks_vecvec(a,b)
  na = length(a)
  nb = length(b)
  na == nb
end

# Matrix vs Vector

function kernel_cache(
  k::FieldBinOp,a::AbstractMatrix,b::AbstractVector)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(k.op,Ta,Tb)
  r = zeros(T,size(a))
  CachedArray(r)
end

function kernel_testitem!(
  c,k::FieldBinOp,a::AbstractMatrix,b::AbstractVector)
  if _valid_checks_matvec(a,b) 
    apply_kernel!(c,k,a,b)
  else
    c.array
  end
end

@inline function apply_kernel!(
  c,k::FieldBinOp,a::AbstractMatrix,b::AbstractVector)
  _field_bin_op_checks_matvec(a,b)
  s = size(a)
  setsize!(c,s)
  np, ni = s
  r = c.array
  for i in 1:ni
    for p in 1:np
      @inbounds r[p,i] = k.op(a[p,i],b[p])
    end
  end
  r
end

function _field_bin_op_checks_matvec(a,b)
  @assert _valid_checks_matvec(a,b) "Binary operation between fields: matrix vs vector size mismatch."
end

function _valid_checks_matvec(a,b)
  na, _ = size(a)
  nb = length(b)
  na == nb
end

# Vector vs Matrix

function kernel_cache(
  k::FieldBinOp,a::AbstractVector,b::AbstractMatrix)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(k.op,Ta,Tb)
  r = zeros(T,size(b))
  CachedArray(r)
end

function kernel_testitem!(
  c,k::FieldBinOp,a::AbstractVector,b::AbstractMatrix)
  if _valid_checks_vecmat(a,b)
    apply_kernel!(c,k,a,b)
  else
    c.array
  end
end

@inline function apply_kernel!(
  c,k::FieldBinOp,a::AbstractVector,b::AbstractMatrix)
  _field_bin_op_checks_vecmat(a,b)
  s = size(b)
  setsize!(c,s)
  np, ni = s
  r = c.array
  for i in 1:ni
    for p in 1:np
      @inbounds r[p,i] = k.op(a[p],b[p,i])
    end
  end
  r
end

function _field_bin_op_checks_vecmat(a,b)
  @assert _valid_checks_vecmat(a,b) "Binary operation between fields: vector vs matrix size mismatch."
end

function _valid_checks_vecmat(a,b)
  nb, _ = size(b)
  na = length(a)
  na == nb
end

# Matrix vs matrix

function kernel_cache(
  k::FieldBinOp,a::AbstractMatrix,b::AbstractMatrix)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(k.op,Ta,Tb)
  np, ni = size(a)
  _, nj = size(b)
  r = zeros(T,(np,ni,nj))
  CachedArray(r)
end

function kernel_testitem!(
  c,k::FieldBinOp,a::AbstractMatrix,b::AbstractMatrix)
  if _valid_checks_matmat(a,b)
    apply_kernel!(c,k,a,b)
  else
    c.array
  end
end

@inline function apply_kernel!(
  c,k::FieldBinOp,a::AbstractMatrix,b::AbstractMatrix)
  _field_bin_op_checks_matmat(a,b)
  np, ni = size(a)
  _, nj = size(b)
  setsize!(c,(np,ni,nj))
  r = c.array
  for j in 1:nj
    for i in 1:ni
      for p in 1:np
        @inbounds r[p,i,j] = k.op(a[p,i],b[p,j])
      end
    end
  end
  r
end

function _field_bin_op_checks_matmat(a,b)
  @assert _valid_checks_matmat(a,b) "Binary operation between fields: matrix vs matrix size mismatch."
end

function _valid_checks_matmat(a,b)
  na, ni = size(a)
  nb, nj = size(b)
  na == nb
end

# Define gradients

function apply_kernel_gradient(k::FieldBinOp,a,b)
  @notimplemented "The gradient of the result of operation $(k.op) is not yet implemented."
end

for op in (:+,:-)
  @eval begin

    function apply_kernel_gradient(k::FieldBinOp{typeof($op)},a,b)
      ga = field_gradient(a)
      gb = field_gradient(b)
      apply_kernel_to_field(k,ga,gb)
    end

    function apply_gradient(k::Valued{FieldBinOp{typeof($op)}},a,b)
      ga = field_array_gradient(a)
      gb = field_array_gradient(b)
      apply(k,ga,gb)
    end

  end
end

