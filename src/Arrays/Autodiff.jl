
function autodiff_array_jacobian(a,i_to_x)
  i_to_xd = lazy_map(DualizeMap(ForwardDiff.jacobian),i_to_x)
  i_to_xdv = VariableArray(i_to_xd)
  i_to_v = a(i_to_xdv)
  i_to_f = FunctorVector(i_to_v)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian),i_to_x)
  i_to_jac = lazy_map(AutoDiffMap(ForwardDiff.jacobian),i_to_f,i_to_x,i_to_cfg)
  i_to_jac
end

function autodiff_array_gradient(a,i_to_x)
  i_to_xd = lazy_map(DualizeMap(ForwardDiff.gradient),i_to_x)
  i_to_xdv = VariableArray(i_to_xd)
  i_to_v = a(i_to_xdv)
  i_to_f = FunctorVector(i_to_v)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient),i_to_x)
  i_to_g = lazy_map(AutoDiffMap(ForwardDiff.gradient),i_to_f,i_to_x,i_to_cfg)
  i_to_g
end

function autodiff_array_hessian(a,i_to_x)
  i_to_xd = lazy_map(DualizeMap(ForwardDiff.hessian),i_to_x)
  i_to_xdv = VariableArray(i_to_xd)
  i_to_v = a(i_to_xdv)
  i_to_f = FunctorVector(i_to_v)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.hessian),i_to_x)
  i_to_g = lazy_map(AutoDiffMap(ForwardDiff.hessian),i_to_f,i_to_x,i_to_cfg)
  i_to_g
end

struct VariableArray{T,N,A} <: AbstractArray{T,N}
  variables::A
  function VariableArray(variables::AbstractArray{T,N}) where {T,N}
    A = typeof(variables)
    new{T,N,A}(variables)
  end
end

Base.size(a::VariableArray) = size(a.variables)
Base.IndexStyle(::Type{VariableArray{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.getindex(a::VariableArray,i::Integer) = a.variables[i]
Base.getindex(a::VariableArray{T,N},i::Vararg{Integer,N}) where {T,N} = a.variables[i...]

struct GetIndexAtX{T}
  x::T
end

function (f::GetIndexAtX)(cache,a::AbstractArray,i)
  getindex!(cache,a,i)
end

function (f::GetIndexAtX)(cache,a::LazyArray,i)
  _cache, index_and_item = cache
  cg, cgi, cf = _cache
  gi = getindex!(cg, a.maps, i)
  fi = map((cj,fj) -> f(cj,fj,i),cf,a.args)
  evaluate!(cgi, gi, fi...)
end

function (f::GetIndexAtX)(cache,a::VariableArray,i)
  f.x
end

struct FunAtIndex{A,B}
  cache::A
  array::B
  i::Int
end

function (f::FunAtIndex)(x)
  k = GetIndexAtX(x)
  k(f.cache,f.array,f.i)
end

struct FunctorVector{T,A} <: AbstractVector{T}
  array::A
  function FunctorVector(array::AbstractArray)
    cache = array_cache(array)
    T = typeof(FunAtIndex(cache,array,1))
    A = typeof(array)
    new{T,A}(array)
  end
end

Base.size(a::FunctorVector) = (length(a.array),)
Base.IndexStyle(::Type{<:FunctorVector}) = IndexLinear()
array_cache(a::FunctorVector) = array_cache(a.array)
function getindex!(cache,a::FunctorVector,i::Integer)
  FunAtIndex(cache,a.array,i)
end
function Base.getindex(a::FunctorVector,i::Integer)
  cache = array_cache(a.array)
  FunAtIndex(cache,a.array,i)
end

struct ConfigMap{F} <: Map
 f::F
end

function return_cache(k::ConfigMap{typeof(ForwardDiff.gradient)},x)
  cfg = ForwardDiff.GradientConfig(nothing,x)
  cfg
end

function return_cache(k::ConfigMap{typeof(ForwardDiff.jacobian)},x)
  cfg = ForwardDiff.JacobianConfig(nothing,x)
  cfg
end

function return_cache(k::ConfigMap{typeof(ForwardDiff.hessian)},x)
  cfg = ForwardDiff.HessianConfig(nothing,x)
  cfg
end

function evaluate!(cfg,k::ConfigMap,x)
  cfg
end

struct DualizeMap{F} <: Map
 f::F
end

function return_cache(k::DualizeMap,x)
  cfg = ConfigMap(k.f)(x)
  xdual = cfg.duals
  xdual
end

function return_cache(k::DualizeMap{typeof(ForwardDiff.hessian)},x)
  cfg = ConfigMap(k.f)(x)
  xdual = cfg.jacobian_config.duals
  xdual
end

function evaluate!(xdual,k::DualizeMap,x)
  @notimplementedif length(xdual) != length(x)
  xdual
end

struct AutoDiffMap{F} <: Map
  f::F
end

function return_cache(k::AutoDiffMap,f,x,cfg)
  k.f(f,x,cfg)
end

function evaluate!(cache,k::AutoDiffMap{typeof(ForwardDiff.jacobian)},f,x,cfg)
  ForwardDiff.jacobian!(cache,f,x,cfg)
  cache
end

function evaluate!(cache,k::AutoDiffMap{typeof(ForwardDiff.gradient)},f,x,cfg)
  ForwardDiff.gradient!(cache,f,x,cfg)
  cache
end

function evaluate!(cache,k::AutoDiffMap{typeof(ForwardDiff.hessian)},f,x,cfg)
  ForwardDiff.hessian!(cache,f,x,cfg)
  cache
end

## OLD
#
#"""
#"""
#function autodiff_array_gradient(a,i_to_x,j_to_i=IdentityVector(length(i_to_x)))
#
#  i_to_xdual = lazy_map(i_to_x) do x
#    cfg = ForwardDiff.GradientConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
#    xdual = cfg.duals
#    xdual
#  end
#
#  j_to_f = to_array_of_functions(a,i_to_xdual,j_to_i)
#  j_to_x = lazy_map(Reindex(i_to_x),j_to_i)
#
#  k = ForwardDiffGradientMap()
#  lazy_map(k,j_to_f,j_to_x)
#
#end
#
#struct ForwardDiffGradientMap <: Map end
#
#function return_cache(k::ForwardDiffGradientMap,f,x)
#  cfg = ForwardDiff.GradientConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
#  r = copy(x)
#  (r, cfg)
#end
#
#@inline function evaluate!(cache,k::ForwardDiffGradientMap,f,x)
#  r, cfg = cache
#  @notimplementedif length(r) != length(x)
#  ForwardDiff.gradient!(r,f,x,cfg)
#  r
#end
#
#"""
#"""
#function autodiff_array_jacobian(a,i_to_x,j_to_i=IdentityVector(length(i_to_x)))
#
#  i_to_xdual = lazy_map(i_to_x) do x
#    cfg = ForwardDiff.JacobianConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
#    xdual = cfg.duals
#    xdual
#  end
#
#  j_to_f = to_array_of_functions(a,i_to_xdual,j_to_i)
#  j_to_x = lazy_map(Reindex(i_to_x),j_to_i)
#
#  k = ForwardDiffJacobianMap()
#  lazy_map(k,j_to_f,j_to_x)
#
#end
#
#struct ForwardDiffJacobianMap <: Map end
#
#function return_cache(k::ForwardDiffJacobianMap,f,x)
#  cfg = ForwardDiff.JacobianConfig(nothing, x, ForwardDiff.Chunk{length(x)}())
#  n = length(x)
#  j = zeros(eltype(x),n,n)
#  (j, cfg)
#end
#
#@inline function evaluate!(cache,k::ForwardDiffJacobianMap,f,x)
#  j, cfg = cache
#  @notimplementedif size(j,1) != length(x)
#  @notimplementedif size(j,2) != length(x)
#  ForwardDiff.jacobian!(j,f,x,cfg)
#  j
#end
#
#"""
#"""
#function autodiff_array_hessian(a,i_to_x,j_to_i=IdentityVector(length(i_to_x)))
#   agrad = i_to_y -> autodiff_array_gradient(a,i_to_y,j_to_i)
#   autodiff_array_jacobian(agrad,i_to_x,j_to_i)
#end
#
#function to_array_of_functions(a,x,ids=IdentityVector(length(x)))
#  k = ArrayOfFunctionsMap(a,x)
#  j = IdentityVector(length(ids))
#  lazy_map(k,j)
#end
#
#struct ArrayOfFunctionsMap{A,X} <: Map
#  a::A
#  x::X
#end
#
#function return_cache(k::ArrayOfFunctionsMap,j)
#  xi = testitem(k.x)
#  l = length(k.x)
#  x = MutableFill(xi,l)
#  ax = k.a(x)
#  axc = array_cache(ax)
#  (ax, x, axc)
#end
#
#@inline function evaluate!(cache,k::ArrayOfFunctionsMap,j)
#  ax, x, axc = cache
#  @inline function f(xj)
#    x.value = xj
#    axj = getindex!(axc,ax,j)
#  end
#  f
#end
#
#mutable struct MutableFill{T} <: AbstractVector{T}
#  value::T
#  length::Int
#end
#
#Base.size(a::MutableFill) = (a.length,)
#
#@inline Base.getindex(a::MutableFill,i::Integer) = a.value
