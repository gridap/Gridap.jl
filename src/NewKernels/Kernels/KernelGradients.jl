# Default gradient and hessian implementation

"""
$(SIGNATURES)

Returns another field that represents the gradient of the given one
"""
function gradient(f)
  FieldGradient(f)
end

struct KernelGradient{F} <: NewKernel
  field::F
  KernelGradient(f) = new{typeof(f)}(f)
end

@inline function evaluate!(cache,f::KernelGradient,x)
  evaluate_gradient!(cache,f.field,x)
end

@inline function return_cache(f::KernelGradient,x)
  return_gradient_cache(f.field,x)
end

@inline function gradient(f::KernelGradient)
  FieldHessian(f.field)
end

struct KernelHessian{F} <: NewKernel
  field::F
  KernelHessian(f) = new{typeof(f)}(f)
end

@inline function evaluate_field!(cache,f::KernelHessian,x)
  evaluate_hessian!(cache,f.field,x)
end

@inline function return_cache(f::KernelHessian,x)
  return_hessian_cache(f.field,x)
end

@inline function gradient(f::KernelHessian)
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
