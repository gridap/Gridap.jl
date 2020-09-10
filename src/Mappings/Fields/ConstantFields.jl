struct ConstantField{T<:Number} <: NewField
  v::T
  function ConstantField{T}(v::T) where {T<:Number}
    new{typeof(v)}(v)
  end
end

ConstantField(v::Number) = ConstantField{typeof(v)}(v)

constant_field(v::Number) =  ConstantField(v)

# Number

function return_cache(f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  c = zeros(typeof(f.v),nx)
  CachedArray(c)
end

function evaluate!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  setsize!(c,(nx,))
  r = c.array
  for i in eachindex(x)
    @inbounds r[i] = f.v
  end
  r
end

# @santiagobadia : It seems wrong to me
function gradient(f::ConstantField)
  T = typeof(f.v)
  E = eltype(T)
  ConstantField(zero(E))
end

# Array

function return_cache(f::ConstantField{<:AbstractArray{<:Number}},x::AbstractArray{<:Point})
  nx = length(x)
  sv = size(f.v)
  s = (nx,sv...)
  c = zeros(eltype(f.v),s)
  CachedArray(c)
end

function evaluate!(c,f::ConstantField{<:AbstractArray{<:Number}},x::AbstractArray{<:Point})
  nx = length(x)
  sv = size(f.v)
  s = (nx,sv...)
  setsize!(c,s)
  cis = CartesianIndices(f.v)
  r = c.array
  for i in eachindex(x)
    for ci in cis
      @inbounds r[i,ci] = f.v[ci]
    end
  end
  r
end

function gradient(f::ConstantField{<:AbstractArray{<:Number}})
  T = eltype(f.v)
  E = eltype(T)
  ConstantField(Fill(zero(E),size(f.v)))
end

using BenchmarkTools
