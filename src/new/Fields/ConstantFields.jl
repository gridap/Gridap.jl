
function field_cache(v::Number,x)
  nx = length(x)
  c = zeros(typeof(v),nx)
  CachedArray(c)
end

function evaluate_field!(c,v::Number,x)
  nx = length(x)
  setsize!(c,(nx,))
  for i in eachindex(x)
    @inbounds c[i] = v
  end
  c
end

function field_gradient(v::Number)
  T = typeof(v)
  E = eltype(T)
  zero(E)
end

function field_cache(v::AbstractArray{<:Number},x)
  nx = length(x)
  sv = size(v)
  s = (nx,sv...)
  c = zeros(eltype(v),s)
  CachedArray(c)
end

function evaluate_field!(c,v::AbstractArray{<:Number},x)
  nx = length(x)
  sv = size(v)
  s = (nx,sv...)
  setsize!(c,s)
  cis = CartesianIndices(v)
  for i in eachindex(x)
    for ci in cis
      @inbounds c[i,ci] = v[ci]
    end
  end
  c
end

function field_gradient(v::AbstractArray{<:Number})
  T = eltype(v)
  E = eltype(T)
  Fill(zero(E),size(v))
end



