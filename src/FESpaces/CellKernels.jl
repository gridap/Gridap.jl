"""
"""
function apply_cellmatvec(op::Function,a,b,c...)
  k = CellMatVecKernel(op)
  apply(k,a,b,c...)
end

"""
"""
function apply_cellmatrix(op::Function,a,b,c...)
  k = CellMatKernel(op)
  apply(k,a,b,c...)
end

"""
"""
function apply_cellvector(op::Function,a,c...)
  k = CellVecKernel(op)
  apply(k,a,c...)
end

# Helpers

struct CellMatVecKernel{F<:Function} <: Kernel
  op::F
end

function kernel_cache(k::CellMatVecKernel,a::AbstractArray,b::AbstractArray,c...)
  @assert ndims(a) == 2
  @assert ndims(b) == 2
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(inner,Ta,Tb)
  m = size(a,2)
  n = size(b,2)
  mat = CachedArray(zeros(T,m,n))
  vec = CachedArray(zeros(T,m))
  (mat, vec)
end

@inline function apply_kernel!(cache,k::CellMatVecKernel,a::AbstractArray,b::AbstractArray,c...)
  cmat, cvec = cache
  m = size(a,2)
  n = size(b,2)
  setsize!(cmat,(m,n))
  setsize!(cvec,(m,))
  mat = cmat.array
  vec = cvec.array
  z = zero(eltype(mat))
  fill!(mat,z)
  fill!(vec,z)
  k.op(mat,vec,a,b,c...)
  (mat,vec)
end

struct CellMatKernel{F<:Function} <: Kernel
  op::F
end

function kernel_cache(k::CellMatKernel,a::AbstractArray,b::AbstractArray,c...)
  @assert ndims(a) == 2
  @assert ndims(b) == 2
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(inner,Ta,Tb)
  m = size(a,2)
  n = size(b,2)
  mat = CachedArray(zeros(T,m,n))
  mat
end

@inline function apply_kernel!(cmat,k::CellMatKernel,a::AbstractArray,b::AbstractArray,c...)
  m = size(a,2)
  n = size(b,2)
  setsize!(cmat,(m,n))
  mat = cmat.array
  z = zero(eltype(mat))
  fill!(mat,z)
  k.op(mat,a,b,c...)
  mat
end

struct CellVecKernel{F<:Function} <: Kernel
  op::F
end

function kernel_cache(k::CellVecKernel,a::AbstractArray,c...)
  @assert ndims(a) == 2
  Ta = eltype(a)
  T = return_type(inner,Ta,Ta)
  m = size(a,2)
  vec = CachedArray(zeros(T,m))
  vec
end

@inline function apply_kernel!(cvec,k::CellVecKernel,a::AbstractArray,c...)
  m = size(a,2)
  setsize!(cvec,(m,))
  vec = cvec.array
  z = zero(eltype(vec))
  fill!(vec,z)
  k.op(vec,a,c...)
  vec
end

