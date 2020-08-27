
"""
    apply_kernel_to_field(k,f...) -> Field

Returns a field obtained by applying the kernel `k` to the 
values of the fields in `f`. That is, the returned field evaluated at
a vector of points `x` provides the value obtained by applying kernel `k` to the
values of the fields `f` at the vector of points `x`. Formally, the resulting field at a 
vector of points
 `x` is defined as

    fx = evaluate_fields(f,x)
    apply_kernel(k,fx...)

In order to be able to call the [`field_gradient`](@ref) function of the
resulting field, one needs to define the gradient operator
associated with the underlying kernel.
This is done by adding a new method to [`apply_kernel_gradient(k,f...)`](@ref) for each kernel type.
"""
@inline function apply_kernel_to_field(k,f...)
  AppliedField(k,f...)
end

"""
    apply_kernel_gradient(k,f...)

Returns a field representing the gradient of the field obtained with

    apply_kernel_to_field(k,f...)
"""
function apply_kernel_gradient(k,f...)
  @abstractmethod
end

# Result of applying a kernel to the value of some fields

struct AppliedField{K,F} <: Field
  k::K
  f::F
  @inline function AppliedField(k,f...)
    new{typeof(k),typeof(f)}(k,f)
  end
end

function field_return_type(f::AppliedField,x)
  Ts = field_return_types(f.f,x)
  kernel_return_type(f.k, testvalues(Ts...)...)
end

function field_cache(f::AppliedField,x)
  cf = field_caches(f.f,x)
  fx = evaluate_fields!(cf,f.f,x)
  ck = kernel_cache(f.k,fx...)
  (ck,cf)
end

@inline function evaluate_field!(cache,f::AppliedField,x)
  ck, cf = cache
  fx = evaluate_fields!(cf,f.f,x)
  apply_kernel!(ck,f.k,fx...)
end

@inline function field_gradient(f::AppliedField)
  apply_kernel_gradient(f.k,f.f...)
end

