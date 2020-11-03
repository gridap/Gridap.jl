
function return_cache(v::Number,x)
  nx = length(x)
  c = zeros(typeof(v),nx)
  CachedArray(c)
end

function evaluate!(c,v::Number,x)
  nx = length(x)
  setsize!(c,(nx,))
  r = c.array
  for i in eachindex(x)
    @inbounds r[i] = v
  end
  r
end

function field_gradient(v::Number)
  T = typeof(v)
  E = eltype(T)
  zero(E)
end

function return_cache(v::AbstractArray{<:Number},x)
  nx = length(x)
  sv = size(v)
  s = (nx,sv...)
  c = zeros(eltype(v),s)
  CachedArray(c)
end

function evaluate!(c,v::AbstractArray{<:Number},x)
  nx = length(x)
  sv = size(v)
  s = (nx,sv...)
  setsize!(c,s)
  cis = CartesianIndices(v)
  r = c.array
  for i in eachindex(x)
    for ci in cis
      @inbounds r[i,ci] = v[ci]
    end
  end
  r
end

function field_gradient(v::AbstractArray{<:Number})
  T = eltype(v)
  E = eltype(T)
  Fill(zero(E),size(v))
end
