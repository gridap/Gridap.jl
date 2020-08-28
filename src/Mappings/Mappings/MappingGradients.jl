# Default gradient and hessian implementation

"""
$(SIGNATURES)

Returns another field that represents the gradient of the given one
"""
function gradient(f)
  FieldGradient(f)
end

struct MappingGradient{F} <: Mapping
  field::F
  MappingGradient(f) = new{typeof(f)}(f)
end

@inline function evaluate!(cache,f::MappingGradient,x)
  evaluate_gradient!(cache,f.field,x)
end

@inline function return_cache(f::MappingGradient,x)
  return_gradient_cache(f.field,x)
end

@inline function gradient(f::MappingGradient)
  FieldHessian(f.field)
end

struct MappingHessian{F} <: Mapping
  field::F
  MappingHessian(f) = new{typeof(f)}(f)
end

@inline function evaluate!(cache,f::MappingHessian,x)
  evaluate_hessian!(cache,f.field,x)
end

@inline function return_cache(f::MappingHessian,x)
  return_hessian_cache(f.field,x)
end

@inline function gradient(f::MappingHessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

"""
$(SIGNATURES)
"""
function evaluate_gradient!(cache,f,x)
  @abstractmethod
end

"""
$(SIGNATURES)
"""
function return_gradient_cache(f,x)
  @abstractmethod
end

"""
$(SIGNATURES)
"""
function evaluate_hessian!(cache,f,x)
  @abstractmethod
end

"""
$(SIGNATURES)
"""
function return_hessian_cache(cache,x)
  @abstractmethod
end
