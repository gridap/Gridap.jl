
"""
"""
function autodiff_array_gradient(a,i_to_x)

  i_to_xdual = apply(i_to_x) do x
    cfg = ForwardDiff.GradientConfig(nothing, x)
    xdual = cfg.duals
    xdual
  end

  i_to_f = to_array_of_functions(a,i_to_xdual)

  k = ForwardDiffGradientKernel()
  apply(k,i_to_f,i_to_x)

end

struct ForwardDiffGradientKernel <: Kernel end

function kernel_cache(k::ForwardDiffGradientKernel,f,x)
  cfg = ForwardDiff.GradientConfig(nothing, x)
  r = copy(x)
  (r, cfg)
end

function apply_kernel!(cache,k::ForwardDiffGradientKernel,f,x)
  r, cfg = cache
  @notimplementedif length(r) != length(x)
  ForwardDiff.gradient!(r,f,x,cfg)
  r
end

"""
"""
function autodiff_array_jacobian(a,i_to_x)

  i_to_xdual = apply(i_to_x) do x
    cfg = ForwardDiff.JacobianConfig(nothing, x)
    xdual = cfg.duals
    xdual
  end

  i_to_f = to_array_of_functions(a,i_to_xdual)

  k = ForwardDiffJacobianKernel()
  apply(k,i_to_f,i_to_x)

end

struct ForwardDiffJacobianKernel <: Kernel end

function kernel_cache(k::ForwardDiffJacobianKernel,f,x)
  cfg = ForwardDiff.JacobianConfig(nothing, x)
  n = length(x)
  j = zeros(eltype(x),n,n)
  (j, cfg)
end

function apply_kernel!(cache,k::ForwardDiffJacobianKernel,f,x)
  j, cfg = cache
  @notimplementedif size(j,1) != length(x)
  @notimplementedif size(j,2) != length(x)
  ForwardDiff.jacobian!(j,f,x,cfg)
  j
end

"""
"""
function autodiff_array_hessian(a,i_to_x)
   agrad = i_to_y -> autodiff_array_gradient(a,i_to_y)
   autodiff_array_jacobian(agrad,i_to_x)
end

function to_array_of_functions(a,x)
  k = ArrayOfFunctionsKernel(a,x)
  i = IdentityVector(length(x))
  apply(k,i)
end

struct ArrayOfFunctionsKernel{A,X} <: Kernel
  a::A
  x::X
end

function kernel_cache(k::ArrayOfFunctionsKernel,i)
  xi = testitem(k.x)
  l = length(k.x)
  x = MutableFill(xi,l)
  ax = k.a(x)
  axc = array_cache(ax)
  (ax, x, axc)
end

function apply_kernel!(cache,k::ArrayOfFunctionsKernel,i)
  ax, x, axc = cache
  @noinline function f(xi)
    x.value = xi
    axi = getindex!(axc,ax,i)
  end
  f
end

mutable struct MutableFill{T} <: AbstractVector{T}
  value::T
  length::Int
end

Base.size(a::MutableFill) = (a.length,)

@inline Base.getindex(a::MutableFill,i::Integer) = a.value

