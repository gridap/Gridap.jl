
"""
    struct MockField{T<:Number} <: Field

For tests.
"""
struct MockField{T<:Number} <: Field
  v::T
end

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

# This is to emulate MonomialBasis that only implements at the level of array of fields.
"""
    struct MockFieldArray{N,A} <: AbstractArray{GenericField{Nothing},N}

For tests.
"""
struct MockFieldArray{N,A} <: AbstractArray{GenericField{Nothing},N}
  values::A
  function MockFieldArray(values::AbstractArray{<:Number})
    N = ndims(values)
    A = typeof(values)
    new{N,A}(values)
  end
end

Base.size(a::MockFieldArray) = size(a.values)
Base.IndexStyle(::Type{<:MockFieldArray}) = IndexLinear()
Base.getindex(a::MockFieldArray,i::Integer) = GenericField(nothing)

function return_cache(f::MockFieldArray,x::Point)
  nothing
end

function evaluate!(c,f::MockFieldArray,x::Point)
  f.values
end

function return_cache(f::MockFieldArray,x::AbstractVector{<:Point})
  T = eltype(f.values)
  CachedArray(zeros(T,length(x),length(f.values)))
end

function evaluate!(c,f::MockFieldArray,x::AbstractVector{<:Point})
  s = (length(x),length(f.values))
  setsize!(c,s)
  r = c.array
  for p in 1:length(x)
    for i in eachindex(f.values)
      r[p,i] = f.values[i]
    end
  end
  r
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

function return_cache(f::FieldGradientArray{1,<:MockFieldArray},x::AbstractVector{<:Point})
  T = return_type(outer,testitem(x),testitem(f.fa.values))
  r = zeros(T,length(x),length(f.fa.values))
  CachedArray(r)
end

function evaluate!(cache,f::FieldGradientArray{1,<:MockFieldArray},x::AbstractVector{<:Point})
  setsize!(cache,(length(x),length(f.fa.values)))
  r = cache.array
  fill!(r,zero(eltype(r)))
  r
end

