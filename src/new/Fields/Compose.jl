
"""
    compose(g::Function,f...)

Returns a new field obtained by composition of function `g` and the fields
`f`. The value of the resulting field at a vector of points `x` is numerically equivalent to

    fx = evaluate_fields(f,x)
    apply_kernel(bcast(g), fx...)

The gradient of the resulting field evaluated at a vector of points `x` is equivalent to

    fx = evaluate_fields(f,x)
    apply_kernel(bcast(gradient(g)), fx...)

Note that it is needed to overload `gradient(::typeof(g))` for the given function `g`
in order to be able to compute the gradient.

As in function [`apply_kernel_to_field`](@ref) if any of the inputs in `f` is a number or an array
instead of a field it will be treated as a "constant field".
"""
function compose(g::Function,f...)
  k = Comp(g)
  apply_kernel_to_field(k,f...)
end

struct Comp{F} <: Kernel
  e::BCasted{F}
  @inline Comp(f::Function) = new{typeof(f)}(BCasted(f))
end

@inline apply_kernel!(cache,k::Comp,x...) = apply_kernel!(cache,k.e,x...)

kernel_cache(k::Comp,x...) = kernel_cache(k.e,x...)

kernel_return_type(k::Comp,x...) = kernel_return_type(k.e,x...)

function apply_kernel_gradient(k::Comp,f...)
  g = gradient(k.e.f)
  compose(g,f...)
end

"""
    compose(g::Function,f::AbstractArray...)

Returns an array of fields numerically equivalent to

    map( (x...)->compose(g,x...), f...)

"""
function compose(g::Function,f::AbstractArray...)
  k = Comp(g)
  apply_to_field_array(k,f...)
end

