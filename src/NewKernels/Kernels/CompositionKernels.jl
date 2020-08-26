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
@inline function composition(k::NewKernel,l::NTuple{N,<:NewKernel}) where N
  CompositionKernel(k,l...)
end

@inline function composition(k::NewKernel,l::NewKernel)
  CompositionKernel(k,l)
end

"""
    apply_kernel_gradient(k,f...)

Returns a field representing the gradient of the field obtained with

    apply_kernel_to_field(k,f...)
"""
# function apply_kernel_gradient(k,f...)
#   @abstractmethod
# end

# @inline apply_kernel_gradient(k::BCasted{typeof(+)},a) = field_gradient(a)

# @inline apply_kernel_gradient(k::BCasted{typeof(-)},a...) = apply_kernel_to_field(k,field_gradients(a...)...)

# @inline apply_kernel_gradient(k::BCasted{typeof(+)},a...) = apply_kernel_to_field(k,field_gradients(a...)...)

# Result of applying a kernel to the value of some fields

struct CompositionKernel{K,L} <: NewKernel
  k::K
  l::L
  @inline function CompositionKernel(k,l...)
    new{typeof(k),typeof(l)}(k,l)
  end
end

function return_type(c::CompositionKernel,x)
  Ts = return_types(c.l,x)
  return_type(c.k, testvalues(Ts...)...)
end

function return_cache(c::CompositionKernel,x)
  cl = caches(c.l,x)
  lx = evaluate!(cl,c.l,x)
  ck = cache(c.k,lx...)
  (ck,cl)
end

@inline function evaluate!(cache,c::CompositionKernel,x)
  ck, cf = cache
  lx = evaluate!(cf,c.l,x)
  evaluate!(ck,c.k,lx...)
end

# @inline function field_gradient(f::CompositionKernel)
  # apply_kernel_gradient(f.k,f.f...)
# end
