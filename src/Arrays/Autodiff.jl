
function autodiff_array_gradient(a,i_to_x)
  dummy_tag = ()->()
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  i_to_ydual = a(i_to_xdual)
  i_to_result = lazy_map(AutoDiffMap(),i_to_cfg,i_to_ydual)
  i_to_result
end

function autodiff_array_jacobian(a,i_to_x)
  dummy_tag = ()->()
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian,dummy_tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  i_to_ydual = a(i_to_xdual)
  i_to_result = lazy_map(AutoDiffMap(),i_to_cfg,i_to_ydual)
  i_to_result
end

function autodiff_array_hessian(a,i_to_x)
  agrad = i_to_y -> autodiff_array_gradient(a,i_to_y)
  autodiff_array_jacobian(agrad,i_to_x)
end

function autodiff_array_gradient(a,i_to_x,j_to_i)
  dummy_tag = ()->()
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_cfg = lazy_map(Reindex(i_to_cfg),j_to_i)
  j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
  j_to_result
end

function autodiff_array_jacobian(a,i_to_x,j_to_i)
  dummy_tag = ()->()
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian,dummy_tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_cfg = lazy_map(Reindex(i_to_cfg),j_to_i)
  j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
  j_to_result
end

function autodiff_array_hessian(a,i_to_x,i_to_j)
  agrad = i_to_y -> autodiff_array_gradient(a,i_to_y,i_to_j)
  autodiff_array_jacobian(agrad,i_to_x,i_to_j)
end

struct ConfigMap{
  F <: Union{typeof(ForwardDiff.gradient),typeof(ForwardDiff.jacobian)},
  T <: Union{<:Function,Nothing}} <: Map

  f::F # ForwardDiff operation
  tag::T # function for config tag name
end

ConfigMap(f) = ConfigMap(f,nothing)

# TODO Prescribing long chunk size can lead to slow compilation times!
function return_cache(k::ConfigMap{typeof(ForwardDiff.gradient)},x)
  cfg = ForwardDiff.GradientConfig(k.tag,x,ForwardDiff.Chunk{length(x)}())
  cfg
end

function return_cache(k::ConfigMap{typeof(ForwardDiff.jacobian)},x)
  cfg = ForwardDiff.JacobianConfig(k.tag,x,ForwardDiff.Chunk{length(x)}())
  cfg
end

function evaluate!(cfg,k::ConfigMap,x)
  cfg
end

struct DualizeMap <: Map end

function evaluate!(cache,::DualizeMap,cfg,x)
  xdual, seeds = cfg.duals, cfg.seeds
  ForwardDiff.seed!(xdual, x, seeds)
  xdual
end

struct AutoDiffMap <: Map end

function return_cache(::AutoDiffMap,cfg::ForwardDiff.GradientConfig,ydual)
  ydual isa Real || throw(ForwardDiff.GRAD_ERROR)
  result = similar(cfg.duals, ForwardDiff.valtype(ydual))
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::ForwardDiff.GradientConfig{T},ydual) where T
  @check ForwardDiff.chunksize(cfg) == length(result)
  result = ForwardDiff.extract_gradient!(T, result, ydual)
  return result
end

function return_cache(::AutoDiffMap,cfg::ForwardDiff.JacobianConfig{T,V,N},ydual) where {T,V,N}
  ydual isa AbstractArray || throw(ForwardDiff.JACOBIAN_ERROR)
  result = similar(ydual, ForwardDiff.valtype(eltype(ydual)), length(ydual), N)
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::ForwardDiff.JacobianConfig{T,V,N},ydual) where {T,V,N}
  @check ForwardDiff.chunksize(cfg) == size(result,2)
  ForwardDiff.extract_jacobian!(T, result, ydual, N)
  ForwardDiff.extract_value!(T, result, ydual)
  return result
end
