
"""
    integrate(f,x,w,j)
"""
function integrate(f,x,w,j)
  fx = evaluate_field(f,x)
  jx = evaluate_field(j,x)
  k = IntKernel()
  lazy_map_kernel(k,fx,w,jx)
end

"""
    integrate(f::AbstractArray,x,w,j)
"""
function integrate(
  f::AbstractArray,x,w,j)
  fx = evaluate_field_array(f,x)
  jx = evaluate_field_array(j,x)
  k = IntKernel()
  lazy_map(k,fx,w,jx)
end

struct IntKernel <: Kernel end

function kernel_cache(k::IntKernel,f::AbstractVector,w,j)
  T = _integrate_rt(f,w,j)
  zero(T)
end

function kernel_testitem!(c,k::IntKernel,f::AbstractVector,w,j)
  if _integrate_valid_sizes(f,w,j)
    lazy_map_kernel!(c,k,f,w,j)
  else
    c
  end
end

@noinline function lazy_map_kernel!(z,k::IntKernel,f::AbstractVector,w,j)
  _integrate_checks(f,w,j)
  r = z
  for p in eachindex(f)
    @inbounds r = r + f[p]*w[p]*meas(j[p])
  end
  r
end

function kernel_cache(k::IntKernel,f::AbstractArray,w,j)
  T = _integrate_rt(f,w,j)
  _, s = _split(size(f)...)
  r = zeros(T,s)
  c = CachedArray(r)
end

function kernel_testitem!(c,k::IntKernel,f::AbstractArray,w,j)
  if _integrate_valid_sizes(f,w,j)
    lazy_map_kernel!(c,k,f,w,j)
  else
    c.array
  end
end

@inline function lazy_map_kernel!(c,k::IntKernel,f::AbstractArray,w,j)
  _integrate_checks(f,w,j)
  np, s = _split(size(f)...)
  cis = CartesianIndices(s)
  setsize!(c,s)
  z = zero(eltype(c))
  r = c.array
  for i in cis
    @inbounds r[i] = z
  end
  for p in 1:np
    @inbounds dV = w[p]*meas(j[p])
    for i in cis
      @inbounds r[i] += f[p,i]*dV
    end
  end
  r
end

function _integrate_rt(f,w,j)
  Tf = eltype(f)
  Tw = eltype(w)
  Tj = eltype(j)
  return_type(*,Tf,Tw,return_type(meas,Tj))
end

function _integrate_checks(f,w,j)
  @assert _integrate_valid_sizes(f,w,j) "integrate: sizes  mismatch."
end

function _integrate_valid_sizes(f,w,j)
  nf, = size(f)
  nw = length(w)
  nj = length(j)
  (nf == nw) && (nw == nj)
end
