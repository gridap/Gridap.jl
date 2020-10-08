
"""
    compose(g::Function,f...)

Returns a new field obtained by composition of function `g` and the fields
`f`. The value of the resulting field at a vector of points `x` is numerically equivalent to

    fx = evaluate_fields(f,x)
    evaluate(bcast(g), fx...)

The gradient of the resulting field evaluated at a vector of points `x` is equivalent to

    fx = evaluate_fields(f,x)
    evaluate(bcast(gradient(g)), fx...)

Note that it is needed to overload `gradient(::typeof(g))` for the given function `g`
in order to be able to compute the gradient.

As in function [`evaluate_to_field`](@ref) if any of the inputs in `f` is a number or an array
instead of a field it will be treated as a "constant field".
"""
function compose(g::Function,f...)
  k = Comp(g)
  evaluate_to_field(k,f...)
end

struct Comp{F} <: Kernel
  e::BCasted{F}
  @inline Comp(f::Function) = new{typeof(f)}(BCasted(f))
end

@inline evaluate!(cache,k::Comp,x...) = evaluate!(cache,k.e,x...)

return_cache(k::Comp,x...) = return_cache(k.e,x...)

return_type(k::Comp,x...) = return_type(k.e,x...)

function evaluate_gradient(k::Comp,f...)
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
  lazy_map_to_field_array(k,f...)
end


"""
    compose_fields(g,f)
"""
function compose_fields(g,f)
  k = CompField(g)
  evaluate_to_field(k,f)
end

"""
    compose(g::Field,f::Field)
"""
function compose(g::Field,f::Field)
  compose_fields(g,f)
end

struct CompField{G} <: Kernel
  g::G
end

return_cache(k::CompField,fx) = field_cache(k.g,fx)

@inline evaluate!(cache,k::CompField,fx) = evaluate_field!(cache,k.g,fx)

return_type(k::CompField,fx) = field_return_type(k.g,fx)

function evaluate_gradient(k::CompField,f)
  g = field_gradient(k.g)
  compose_fields(g,f)
end

"""
    compose_field_arrays(g,f)
"""
function compose_field_arrays(g,f)
  k = CompFieldArray()
  lazy_map(k,g,f)
end

"""
    compose(g::AbstractArray{<:Field},f::AbstractArray{<:Field})
"""
function compose(g::AbstractArray{<:Field},f::AbstractArray{<:Field})
  compose_field_arrays(g,f)
end

struct CompFieldArray <: Kernel end

return_cache(k::CompFieldArray,gi,fi) = nothing

@inline evaluate!(cache,k::CompFieldArray,gi,fi) = compose_fields(gi,fi)

function kernel_evaluate(k::CompFieldArray,x,g,f)
  fx = evaluate_field_array(f,x)
  gx = evaluate_field_array(g,fx)
  gx
end

function lazy_map_gradient(k::CompFieldArray,g,f)
  ∇g = field_array_gradient(g)
  compose_field_arrays(∇g,f)
end
