
"""
    varinner(a,b)
"""
function varinner(a,b)
  k = Varinner()
  apply_kernel_to_field(k,a,b)
end

"""
    varinner(a::AbstractArray,b::AbstractArray)
"""
function varinner(a::AbstractArray,b::AbstractArray)
  k = Varinner()
  apply_to_field_array(k,a,b)
end

struct Varinner <: Kernel end

function apply_kernel_gradient(::Varinner,a,b)
  @unreachable "The gradient of the result of varinner does not make sense."
end

# Vector vs Vector

function kernel_cache(
  k::Varinner,a::AbstractVector,b::AbstractVector)
  _varinner_checks_vecvec(a,b)
  na = length(a)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(inner,Ta,Tb)
  r = zeros(T,na)
  CachedArray(r)
end

@inline function apply_kernel!(
  r,k::Varinner,a::AbstractVector,b::AbstractVector)
  _varinner_checks_vecvec(a,b)
  na = length(a)
  setsize!(r,(na,))
  for p in eachindex(a)
    @inbounds r[p] = inner(a[p],b[p])
  end
  r
end

function _varinner_checks_vecvec(a,b)
  na = length(a)
  nb = length(b)
  @assert na == nb "varinner: vector vs vector size mismatch."
end

# Matrix vs Vector

function kernel_cache(
  k::Varinner,a::AbstractMatrix,b::AbstractVector)
  _varinner_checks_matvec(a,b)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(inner,Ta,Tb)
  r = zeros(T,size(a))
  CachedArray(r)
end

@inline function apply_kernel!(
  r,k::Varinner,a::AbstractMatrix,b::AbstractVector)
  _varinner_checks_matvec(a,b)
  s = size(a)
  setsize!(r,s)
  np, ni = s
  for i in 1:ni
    for p in 1:np
      @inbounds r[p,i] = inner(a[p,i],b[p])
    end
  end
  r
end

function _varinner_checks_matvec(a,b)
  na, _ = size(a)
  nb = length(b)
  @assert na == nb "varinner: matrix vs vector size mismatch."
end

# Matrix vs matrix

function kernel_cache(
  k::Varinner,a::AbstractMatrix,b::AbstractMatrix)
  _varinner_checks_matmat(a,b)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(inner,Ta,Tb)
  np, ni = size(a)
  _, nj = size(b)
  r = zeros(T,(np,ni,nj))
  CachedArray(r)
end

@inline function apply_kernel!(
  r,k::Varinner,a::AbstractMatrix,b::AbstractMatrix)
  _varinner_checks_matmat(a,b)
  np, ni = size(a)
  _, nj = size(b)
  setsize!(r,(np,ni,nj))
  for j in 1:nj
    for i in 1:ni
      for p in 1:np
        @inbounds r[p,i,j] = inner(a[p,i],b[p,j])
      end
    end
  end
  r
end

function _varinner_checks_matmat(a,b)
  na, ni = size(a)
  nb, nj = size(b)
  @assert na == nb "varinner: matrix vs matrix size mismatch."
end

