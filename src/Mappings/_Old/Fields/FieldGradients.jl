gradient(f::NewField) = GenericField(Gradient(f))

struct Gradient{F} <: Mapping
  object::F
end

@inline evaluate!(cache,f::Gradient,x) = evaluate_gradient!(cache,f.object,x)

@inline return_cache(f::Gradient,x) = return_gradient_cache(f.object,x)

struct Hessian{F}
  object::F
end

@inline evaluate!(cache,f::Hessian,x) = evaluate_hessian!(cache,f.object,x)

@inline return_cache(f::Hessian,x) = return_hessian_cache(f.object,x)

@inline function gradient(f::Hessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

# Array of gradients

gradient(a::FieldArray) = FieldGradientArray(Gradient(a))

Base.size(a::Gradient) = size(a.object)
Base.getindex(a::FieldArray,i...) = GenericField(Gradient(a.object[i...]))
Base.IndexStyle(::Type{<:Gradient{F}) where F = IndexStyle(F)
Base.ndims(a::Gradient{F}) where F = ndims(F)
Base.eltype(a::Gradient{F}) where F = GenericField{Gradient{eltype(F)}}

function gradient(f::FieldGradientArray)
  FieldHessianArray(f.i_to_f)
end
