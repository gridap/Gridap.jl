
"""
"""
function autodiff_array_gradient(a,i_to_x,j_to_i=IdentityVector(length(i_to_x)))

  i_to_xdual = lazy_map(i_to_x) do x
    cfg = ForwardDiff.GradientConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
    xdual = cfg.duals
    xdual
  end

  j_to_f = to_array_of_functions(a,i_to_xdual,j_to_i)
  j_to_x = reindex(i_to_x,j_to_i)

  k = ForwardDiffGradientMap()
  lazy_map(k,j_to_f,j_to_x)

end

struct ForwardDiffGradientMap <: Map end

function return_cache(k::ForwardDiffGradientMap,f,x)
  cfg = ForwardDiff.GradientConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
  r = copy(x)
  (r, cfg)
end

@inline function evaluate!(cache,k::ForwardDiffGradientMap,f,x)
  r, cfg = cache
  @notimplementedif length(r) != length(x)
  ForwardDiff.gradient!(r,f,x,cfg)
  r
end

"""
"""
function autodiff_array_jacobian(a,i_to_x,j_to_i=IdentityVector(length(i_to_x)))

  i_to_xdual = lazy_map(i_to_x) do x
    cfg = ForwardDiff.JacobianConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
    xdual = cfg.duals
    xdual
  end

  j_to_f = to_array_of_functions(a,i_to_xdual,j_to_i)
  j_to_x = reindex(i_to_x,j_to_i)

  k = ForwardDiffJacobianMap()
  lazy_map(k,j_to_f,j_to_x)

end

struct ForwardDiffJacobianMap <: Map end

function return_cache(k::ForwardDiffJacobianMap,f,x)
  cfg = ForwardDiff.JacobianConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
  n = length(x)
  j = zeros(eltype(x),n,n)
  (j, cfg)
end

@inline function evaluate!(cache,k::ForwardDiffJacobianMap,f,x)
  j, cfg = cache
  @notimplementedif size(j,1) != length(x)
  @notimplementedif size(j,2) != length(x)
  ForwardDiff.jacobian!(j,f,x,cfg)
  j
end

"""
"""
function autodiff_array_hessian(a,i_to_x,j_to_i=IdentityVector(length(i_to_x)))
   agrad = i_to_y -> autodiff_array_gradient(a,i_to_y,j_to_i)
   autodiff_array_jacobian(agrad,i_to_x,j_to_i)
end

function to_array_of_functions(a,x,ids=IdentityVector(length(x)))
  k = ArrayOfFunctionsMap(a,x)
  j = IdentityVector(length(ids))
  lazy_map(k,j)
end

struct ArrayOfFunctionsMap{A,X} <: Map
  a::A
  x::X
end

function return_cache(k::ArrayOfFunctionsMap,j)
  xi = testitem(k.x)
  l = length(k.x)
  x = MutableFill(xi,l)
  ax = k.a(x)
  axc = array_cache(ax)
  (ax, x, axc)
end

@inline function evaluate!(cache,k::ArrayOfFunctionsMap,j)
  ax, x, axc = cache
  @inline function f(xj)
    x.value = xj
    axj = getindex!(axc,ax,j)
  end
  f
end

mutable struct MutableFill{T} <: AbstractVector{T}
  value::T
  length::Int
end

Base.size(a::MutableFill) = (a.length,)

@inline Base.getindex(a::MutableFill,i::Integer) = a.value
