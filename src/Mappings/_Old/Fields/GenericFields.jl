struct GenericField{F} <: NewField
  object::F
end

# Convert to a field if needed
Field(f) = GenericField(f)
GenericField(f::NewField) = f
return_type(::Type{<:GenericField},::Type{T}) where T<:Field = T
return_type(::Type{<:GenericField},::Type{T}) where T = GenericField{T}


return_cache(a::GenericField,x) = return_cache(a.object,x)
evaluate!(cache,a::GenericField,x) = evaluate!(cache,a.object,x)
gradient(a::GenericField) = GenericField(Gradient(a.object))

# Wrap anything that behaves like an array of fields
# field_array needs to implement:
# size, axes, getindex, eltype, ndims, IndexStyle
# Here we assume that the entire collection field_array
# implements efficient `evaluate!` (a default implementation is given above)
struct GenericFieldArray{A,F,N} <: AbstractArray{F,N}
  field_array::A
  function GenericFieldArray(a::AbstractArray{T,N}) where {T<:NewField,N}
    # T = eltype(a)
    # N = ndims(a)
    # F = GenericField{T}
    # new{F,N,typeof(a)}(a)
    new{typeof(a)T,N}(a)
  end
end

# Convert to Field array when needed
FieldArray(a) = GenericFieldArray(a)
GenericFieldArray(a::AbstractArray{<:Field}) = a

return_cache(a::GenericFieldArray,x) = return_cache(a.field_array,x)
evaluate!(cache,a::GenericFieldArray,x) = evaluate!(cache,a.field_array,x)
gradient(a::GenericFieldArray) = GenericFieldArray(Gradient(a.field_array))

Base.size(a::GenericFieldArray) = size(a.field_array)
Base.axes(a::GenericFieldArray) = axes(a.field_array)
Base.getindex(a::GenericFieldArray,i::Integer...) = GenericField(a.field_array[i...])
Base.IndexStyle(::Type{<:GenericFieldArray{T,N,A}}) where {T,N,A} = IndexStyle(A)
