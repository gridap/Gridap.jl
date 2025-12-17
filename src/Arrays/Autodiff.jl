
function autodiff_array_gradient(a,i_to_x)
  tag = x->ForwardDiff.gradient(a, x)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  i_to_ydual = a(i_to_xdual)
  i_to_result = lazy_map(AutoDiffMap(),i_to_cfg,i_to_ydual)
  i_to_result
end

function autodiff_array_jacobian(a,i_to_x)
  tag = x->Forwarddiff.jacobian(a, x)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian,tag),i_to_x)
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
  tag = x->ForwardDiff.gradient(a, x)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_cfg = autodiff_array_reindex(i_to_cfg,j_to_i)
  j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
  j_to_result
end

function autodiff_array_jacobian(a,i_to_x,j_to_i)
  tag = x->Forwarddiff.jacobian(a, x)
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian,tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_cfg = autodiff_array_reindex(i_to_cfg,j_to_i)
  j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
  j_to_result
end

function autodiff_array_hessian(a,i_to_x,j_to_i)
  agrad = i_to_y -> autodiff_array_gradient(a,i_to_y,j_to_i)
  autodiff_array_jacobian(agrad,i_to_x,j_to_i)
end

function autodiff_array_reindex(i_to_val, j_to_i)
  n_neg = count(j -> j < 0, j_to_i)
  if iszero(n_neg)
    j_to_val = lazy_map(Reindex(i_to_val),j_to_i)
  else
    neg_to_val = Fill(testitem(i_to_val),n_neg)
    j_to_val = lazy_map(PosNegReindex(i_to_val,neg_to_val),j_to_i)
  end
  return j_to_val
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
  result = CachedArray(similar(cfg.duals, ForwardDiff.valtype(ydual)))
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::ForwardDiff.GradientConfig{T},ydual) where T
  @check ForwardDiff.chunksize(cfg) == length(result)
  setsize!(result, (ForwardDiff.chunksize(cfg),))
  result = ForwardDiff.extract_gradient!(T, result, ydual)
  return result
end

function return_cache(::AutoDiffMap,cfg::ForwardDiff.JacobianConfig{T,V,N},ydual) where {T,V,N}
  ydual isa AbstractArray || throw(ForwardDiff.JACOBIAN_ERROR)
  result = CachedArray(similar(ydual, ForwardDiff.valtype(eltype(ydual)), length(ydual), N))
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::ForwardDiff.JacobianConfig{T,V,N},ydual) where {T,V,N}
  @check ForwardDiff.chunksize(cfg) == size(result,2)
  if !isempty(ydual)  # TODO: Temporary fix, sometimes ydual.touched is incorrectly true for the
                      #       case of SkeletonTriangulation + MultiField on different triangulations.
    setsize!(result, (length(ydual),N))
    ForwardDiff.extract_jacobian!(T, result, ydual, N)
    ForwardDiff.extract_value!(T, result, ydual)
  end
  return result
end

# Autodiff - Skeleton + SingleField

function return_cache(k::AutoDiffMap,cfg::ForwardDiff.JacobianConfig,ydual::VectorBlock)
  i = findfirst(ydual.touched)
  yi = ydual.array[i]
  ci = return_cache(k,cfg,yi)
  ri = evaluate!(ci,k,cfg,yi)
  cache = Vector{typeof(ci)}(undef,length(ydual.array))
  array = Vector{typeof(ri)}(undef,length(ydual.array))
  for i in eachindex(ydual.array)
    if ydual.touched[i]
      cache[i] = return_cache(k,cfg,ydual.array[i])
    end
  end
  result = ArrayBlock(array,ydual.touched)
  return result, cache
end

function evaluate!(cache,k::AutoDiffMap,cfg::ForwardDiff.JacobianConfig,ydual::VectorBlock)
  r, c = cache
  for i in eachindex(ydual.array)
    if ydual.touched[i]
      r.array[i] = evaluate!(c[i],k,cfg,ydual.array[i])
    end
  end
  return r
end

# Autodiff - MultiField

struct BlockConfig{C,T,V,N,D,O} <: ForwardDiff.AbstractConfig{N}
  seeds::NTuple{N,ForwardDiff.Partials{N,V}}
  duals::D
  offsets::O
end

function BlockConfig(
  op::Union{typeof(ForwardDiff.gradient),typeof(ForwardDiff.jacobian)},
  f::F,
  x::Union{VectorBlock{<:AbstractArray{V}},VectorBlock{<:VectorBlock{<:AbstractArray{V}}}},
  ::T = ForwardDiff.Tag(f,V)
) where {F,V,T}
  offsets, N = block_offsets(x, 0)
  seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})
  duals = similar(x, ForwardDiff.Dual{T,V,N})
  BlockConfig{typeof(op),T,V,N,typeof(duals),typeof(offsets)}(seeds,duals,offsets)
end

@inline block_offsets(x::Vector, offset) = offset, offset + length(x)

function block_offsets(x::VectorBlock, offset)
  offsets = ()
  for i in eachindex(x.touched)
    if x.touched[i]
      @inbounds offsets_i, offset = block_offsets(x.array[i], offset)
    else
      offsets_i = -1
    end
    offsets = (offsets...,offsets_i)
  end
  return offsets, offset
end

for F in (ForwardDiff.gradient,ForwardDiff.jacobian)
  @eval begin
    function return_cache(k::ConfigMap{typeof($F)},x::VectorBlock)
      return BlockConfig($F,k.tag,x)
    end
  end
end

function evaluate!(cache,k::DualizeMap,cfg::BlockConfig,x)
  xdual, seeds, offsets = cfg.duals, cfg.seeds, cfg.offsets
  seed_block!(xdual, x, seeds, offsets)
  return xdual
end

function return_cache(::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.gradient),T},ydual) where T
  ydual isa Real || throw(ForwardDiff.GRAD_ERROR)
  result = CachedArray(similar(cfg.duals, ForwardDiff.valtype(ydual)))
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.gradient),T},ydual) where T
  _setsize!(result, cfg.duals)
  extract_gradient_block!(T, result, ydual, cfg.offsets)
  return result
end

function _setsize!(result::VectorBlock,duals::VectorBlock{<:Vector})
  ni = size(result.array,1)
  for i in 1:ni
    if result.touched[i]
      setsize!(result[i], (length(duals[i]),))
    end
  end
end

function return_cache(::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.jacobian),T},ydual) where T
  ydual isa VectorBlock || throw(ForwardDiff.JACOBIAN_ERROR)
  result = _alloc_jacobian(ydual,cfg.duals)
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.jacobian),T},ydual) where T
  _setsize!(result,ydual)
  extract_jacobian_block!(T, result, ydual, cfg.offsets)
  return result
end

function _alloc_jacobian(ydual::AbstractVector,xdual::AbstractVector)
  T = ForwardDiff.valtype(eltype(ydual))
  CachedArray(zeros(T,length(ydual),length(xdual)))
end

function _setsize!(result::MatrixBlock,ydual::VectorBlock{<:AbstractVector})
  ni,nj = size(result)
  for i in 1:ni
    for j in 1:nj
      if result.touched[i,j]
        setsize!(result[i,j], (length(ydual[i]),length(ydual[j])))
      end
    end
  end
end

# Skeleton + Multifield: The VectorBlock corresponds to +/-
function _alloc_jacobian(ydual::VectorBlock,xdual::AbstractVector)
  i = findfirst(ydual.touched)
  ai = _alloc_jacobian(ydual.array[i],xdual)
  ni = size(ydual.array,1)
  array = Vector{typeof(ai)}(undef,ni)
  for i in 1:ni
    if ydual.touched[i]
      array[i] = _alloc_jacobian(ydual.array[i],xdual)
    end
  end
  ArrayBlock(array,ydual.touched)
end

function _alloc_jacobian(ydual::VectorBlock,xdual::VectorBlock)
  i = findfirst(ydual.touched)
  j = findfirst(xdual.touched)
  ai = _alloc_jacobian(ydual.array[i],xdual.array[j])

  ni, nj = size(ydual.array,1), size(xdual.array,1)
  array = Matrix{typeof(ai)}(undef,ni,nj)
  touched = fill(false,ni,nj)
  for i in 1:ni
    for j in 1:nj
      if ydual.touched[i] && xdual.touched[j]
        array[i,j]   = _alloc_jacobian(ydual.array[i],xdual.array[j])
        touched[i,j] = true
      end
    end
  end
  ArrayBlock(array,touched)
end

function _setsize!(result::MatrixBlock{<:VectorBlock},ydual::VectorBlock{<:VectorBlock})
  ni,nj = size(result)
  for i in 1:ni
    for j in 1:nj
      if result.touched[i,j]
        setsize!(result[i,j][1], (length(ydual[i][1]),length(ydual[j][1])))
        setsize!(result[i,j][2], (length(ydual[i][2]),length(ydual[j][2])))
      end
    end
  end
end

function seed_block!(
  duals::VectorBlock, x::VectorBlock, seeds::NTuple{N,ForwardDiff.Partials{N}}, offsets
) where N
  for i in eachindex(duals.touched)
    if duals.touched[i]
      @check x.touched[i]
      @inbounds seed_block!(duals.array[i], x.array[i], seeds, offsets[i])
    end
  end
  return duals
end

function seed_block!(
  duals::AbstractArray{ForwardDiff.Dual{T,V,N}}, x, seeds::NTuple{N,ForwardDiff.Partials{N,V}}, offset
) where {T,V,N}
  for j in eachindex(duals)
    @inbounds duals[j] = ForwardDiff.Dual{T,V,N}(x[j], seeds[j+offset])
  end
  return duals
end

function extract_gradient_block!(::Type{T}, result::VectorBlock, dual, offsets) where T
  for i in eachindex(result.touched)
    if result.touched[i]
      @inbounds extract_gradient_block!(T, result.array[i], dual, offsets[i])
    end
  end
  return result
end

function extract_gradient_block!(::Type{T}, result::AbstractArray, dual::ForwardDiff.Dual, offset) where {T}
  for j in eachindex(result)
    @inbounds result[j] = ForwardDiff.partials(T,dual,j+offset)
  end
  return result
end

function extract_gradient_block!(::Type{T}, result::AbstractArray, dual::Real, offset) where {T}
  fill!(result,zero(dual))
  return result
end

function extract_jacobian_block!(::Type{T}, result::MatrixBlock, dual::VectorBlock, offsets) where T
  for i in axes(result.touched,1)
    for j in axes(result.touched,2)
      if result.touched[i,j]
        @inbounds extract_jacobian_block!(T, result.array[i,j], dual.array[i], offsets[j])
      end
    end
  end
  return result
end

# Skeleton + Multifield: The VectorBlocks correspond to +/-
function extract_jacobian_block!(::Type{T}, result::VectorBlock, dual::VectorBlock, offset) where T
  for i in axes(result.touched,1)
    if result.touched[i]
      @check dual.touched[i]
      @inbounds extract_jacobian_block!(T, result.array[i], dual.array[i], offset)
    end
  end
  return result
end

function extract_jacobian_block!(::Type{T}, result::AbstractArray, dual::AbstractArray{<:ForwardDiff.Dual}, offset) where {T}
  for k in axes(result,1)
    for l in axes(result,2)
      @inbounds result[k,l] = ForwardDiff.partials(T,dual[k],l+offset)
    end
  end
  return result
end

function extract_jacobian_block!(::Type{T}, result::AbstractArray, dual::AbstractArray{<:Real}, offset) where {T}
  fill!(result,zero(eltype(dual)))
  return result
end
