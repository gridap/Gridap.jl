
"""
    evaluate_to_field(k,f...) -> Field

Returns a field obtained by applying the kernel `k` to the
values of the fields in `f`. That is, the returned field evaluated at
a vector of points `x` provides the value obtained by applying kernel `k` to the
values of the fields `f` at the vector of points `x`. Formally, the resulting field at a
vector of points
 `x` is defined as

    fx = evaluate_fields(f,x)
    evaluate(k,fx...)

In order to be able to call the [`field_gradient`](@ref) function of the
resulting field, one needs to define the gradient operator
associated with the underlying kernel.
This is done by adding a new method to [`evaluate_gradient(k,f...)`](@ref) for each kernel type.
"""
@inline function evaluate_to_field(k,f...)
  AppliedField(k,f...)
end

"""
    evaluate_gradient(k,f...)

Returns a field representing the gradient of the field obtained with

    evaluate_to_field(k,f...)
"""
function evaluate_gradient(k,f...)
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
  return_type(f.k,map(testvalue,Ts)...)
end

function return_cache(f::AppliedField,x)
  cf = return_caches(f.f,x)
  fx = evaluate_fields!(cf,f.f,x)
  ck = return_cache(f.k,fx...)
  (ck,cf)
end

@inline function evaluate!(cache,f::AppliedField,x)
  ck, cf = cache
  fx = evaluate_fields!(cf,f.f,x)
  evaluate!(ck,f.k,fx...)
end

@inline function field_gradient(f::AppliedField)
  evaluate_gradient(f.k,f.f...)
end
