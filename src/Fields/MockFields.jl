
struct MockField{T<:Number} <: Field
  v::T
end

#MockField(D::Int,v::Number) = MockField{D}(v)

#mock_field(D::Int,v::Number) = MockField{D}(v)

function evaluate!(c,f::MockField,x::Point)
  f.v
end

function evaluate!(cache,f::FieldGradient{1,<:MockField},x::Point)
  zero(outer(x,f.object.v))
end

function evaluate!(cache,f::FieldGradient{2,<:MockField},x::Point)
  zero(outer(x,outer(x,f.object.v)))
end

testvalue(::Type{MockField{T}}) where T = MockField(zero(T))

struct MockFieldArray{T,N,A} <: AbstractArray{MockField{T},N}
  values::A
  function MockFieldArray(values::AbstractArray{<:Number})
    T = eltype(values)
    N = ndims(values)
    A = typeof(values)
    new{T,N,A}(values)
  end
end

Base.size(a::MockFieldArray) = size(a.values)
Base.IndexStyle(::Type{<:MockFieldArray}) = IndexLinear()
@inline Base.getindex(a::MockFieldArray,i::Integer) = MockField(a.values[i])

function evaluate!(c,f::MockFieldArray,x::Point)
  f.values
end

function return_cache(f::FieldGradientArray{1,<:MockFieldArray},x::Point)
  T = return_type(outer,x,testitem(f.fa.values))
  r = zeros(T,size(f.fa.values))
  CachedArray(r)
end

function evaluate!(cache,f::FieldGradientArray{1,<:MockFieldArray},x::Point)
  if size(cache) != size(f.fa.values)
    setsize!(cache,size(f.fa.values))
    r = cache.array
    fill!(r,zero(eltype(r)))
  end
  cache.array
end

