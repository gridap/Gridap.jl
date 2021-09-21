
function autodiff_array_gradient(a,i_to_x)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(ForwardDiff.gradient),i_to_x)
  i_to_ydual = a(i_to_xdual)
  i_to_result = lazy_map(AutoDiffMap(ForwardDiff.gradient),i_to_ydual,i_to_x,i_to_cfg)
  i_to_result
end

function autodiff_array_jacobian(a,i_to_x)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(ForwardDiff.jacobian),i_to_x)
  i_to_ydual = a(i_to_xdual)
  i_to_result = lazy_map(AutoDiffMap(ForwardDiff.jacobian),i_to_ydual,i_to_x,i_to_cfg)
  i_to_result
end

function autodiff_array_hessian(a,i_to_x)
  agrad = i_to_y -> autodiff_array_gradient(a,i_to_y)
  autodiff_array_jacobian(agrad,i_to_x)
end

function autodiff_array_gradient(a,i_to_x,j_to_i)
  i_to_xdual = lazy_map(DualizeMap(ForwardDiff.gradient),i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_x = lazy_map(Reindex(i_to_x),j_to_i)
  j_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient),j_to_x)
  j_to_result = lazy_map(AutoDiffMap(ForwardDiff.gradient),j_to_ydual,j_to_x,j_to_cfg)
  j_to_result
end

function autodiff_array_jacobian(a,i_to_x,j_to_i)
  i_to_xdual = lazy_map(DualizeMap(ForwardDiff.jacobian),i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_x = lazy_map(Reindex(i_to_x),j_to_i)
  j_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian),j_to_x)
  j_to_result = lazy_map(AutoDiffMap(ForwardDiff.jacobian),j_to_ydual,j_to_x,j_to_cfg)
  j_to_result
end

function autodiff_array_hessian(a,i_to_x,i_to_j)
  agrad = i_to_y -> autodiff_array_gradient(a,i_to_y,i_to_j)
  autodiff_array_jacobian(agrad,i_to_x,i_to_j)
end

struct ConfigMap{F} <: Map
 f::F
end

# TODO Prescribing long chunk size can lead to slow compilation times!
function return_cache(k::ConfigMap{typeof(ForwardDiff.gradient)},x)
  cfg = ForwardDiff.GradientConfig(nothing,x,ForwardDiff.Chunk{length(x)}())
  cfg
end

function return_cache(k::ConfigMap{typeof(ForwardDiff.jacobian)},x)
  cfg = ForwardDiff.JacobianConfig(nothing,x,ForwardDiff.Chunk{length(x)}())
  cfg
end

function evaluate!(cfg,k::ConfigMap,x)
  cfg
end

struct DualizeMap{F} <: Map
 f::F
end

function return_cache(k::DualizeMap,x)
  return_cache(ConfigMap(k.f),x)
end

function evaluate!(cfg,k::DualizeMap,x)
  xdual = cfg.duals
  ForwardDiff.seed!(xdual, x, cfg.seeds)
  xdual
end

struct AutoDiffMap{F} <: Map
  f::F
end

function return_cache(k::AutoDiffMap,ydual,x,cfg::ForwardDiff.GradientConfig{T}) where T
  ydual isa Real || throw(ForwardDiff.GRAD_ERROR)
  result = similar(x, ForwardDiff.valtype(ydual))
  result
end

function evaluate!(result,k::AutoDiffMap,ydual,x,cfg::ForwardDiff.GradientConfig{T}) where T
  @notimplementedif ForwardDiff.chunksize(cfg) != length(x)
  @notimplementedif length(result) != length(x)
  result = ForwardDiff.extract_gradient!(T, result, ydual)
  return result
end

function return_cache(k::AutoDiffMap,ydual,x,cfg::ForwardDiff.JacobianConfig{T,V,N}) where {F,T,V,N}
  ydual isa AbstractArray || throw(ForwardDiff.JACOBIAN_ERROR)
  result = similar(ydual, ForwardDiff.valtype(eltype(ydual)), length(ydual), N)
  result
end

function evaluate!(result,k::AutoDiffMap,ydual,x,cfg::ForwardDiff.JacobianConfig{T,V,N}) where {F,T,V,N}
  @notimplementedif ForwardDiff.chunksize(cfg) != length(x)
  @notimplementedif size(result,2) != length(x)
  ForwardDiff.extract_jacobian!(T, result, ydual, N)
  ForwardDiff.extract_value!(T, result, ydual)
  return result
end

