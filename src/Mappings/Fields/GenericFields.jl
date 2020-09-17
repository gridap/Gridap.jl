struct GenericField{F} <: NewField
  field::F
end

GenericField(f::NewField) = f


return_cache(a::GenericField,x) = return_cache(a.field,x)
evaluate!(cache,a::GenericField,x) = evaluate!(cache,a.field,x)
gradient(a::GenericField) = GenericField(Gradient(a.field))

struct GenericFieldArray{F,N,A} <: AbstractArray{F,N}
  field_array::A
  function GenericFieldArray(a::AbstractArray{T,N}) where {T<:NewField,N}
    # T = eltype(a)
    # N = ndims(a)
    # F = GenericField{T}
    # new{F,N,typeof(a)}(a)
    new{T,N,typeof(a)}(a)
  end
end

return_cache(a::GenericFieldArray,x) = return_cache(a.field_array,x)
evaluate!(cache,a::GenericFieldArray,x) = evaluate!(cache,a.field_array,x)
gradient(a::GenericFieldArray) = GenericFieldArray(Gradient(a.field_array))

Base.size(a::GenericFieldArray) = size(a.field_array)
Base.getindex(a::GenericFieldArray,i::Integer...) = GenericField(a.field_array[i...])
Base.IndexStyle(::Type{<:GenericFieldArray{T,N,A}}) where {T,N,A} = IndexStyle(A)
