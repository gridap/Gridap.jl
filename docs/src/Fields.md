```@meta
CurrentModule = Gridap.Fields
```
# Gridap.Fields

```@docs
Fields
```

## Interface

```@docs
Field
Point
evaluate_field!(cache,f,x)
field_cache(f,x)
field_gradient(f)
field_return_type(f,x)
evaluate_gradient!(cache,f,x)
gradient_cache(f,x)
evaluate_hessian!(cache,f,x)
hessian_cache(f,x)
test_field
```
## Helper functions using fields

```@docs
evaluate_field(f,x)
evaluate(f::Field,x)
evaluate!(cache,f::Field,x)
gradient(f::Field)
∇
gradient_type
```

## Working with several fields at once

```@docs
field_return_types(f::Tuple,x)
field_caches(f::Tuple,x)
evaluate_fields(f::Tuple,x)
evaluate_fields!(cf::Tuple,f::Tuple,x)
field_gradients(a,b...)
gradient_all(a,b...)
evaluate_all(f::Tuple,x)
```

## Applying kernels to fields

```@docs
lazy_map_kernel_to_field(k,f...)
lazy_map_kernel_gradient(k,f...)
```

## Working with arrays of fields

```@docs
evaluate_field_array(a::AbstractArray,x::AbstractArray)
evaluate(::AbstractArray{<:Field},::AbstractArray)
field_array_gradient(a::AbstractArray)
gradient(::AbstractArray{<:Field})
field_array_cache(a::AbstractArray,x::AbstractArray)
test_array_of_fields
```

## Working with several arrays of fields at once

```@docs
evaluate_field_arrays(f::Tuple,x::AbstractArray)
field_array_gradients(f...)
```
## Applying kernels to arrays of fields

```@docs
lazy_map_to_field_array(k,f::AbstractArray...)
kernel_evaluate(k,x,f...)
lazy_map_gradient(k,f...)
```

## Operations on fields and arrays of fields


```@docs
field_operation(op::Function,a,b)
field_array_operation(op::Function,a,b)
compose_fields(g,f)
compose_field_arrays(g,f)
compose(g::Field,f::Field)
compose(g::AbstractArray{<:Field},f::AbstractArray{<:Field})
compose(g::Function,f...)
compose(g::Function,f::AbstractArray...)
lincomb(a::Field,b::AbstractVector)
lincomb(a::AbstractArray,b::AbstractArray)
lazy_map_lincomb(ax,b)
attachmap(f,phi)
attachmap(f::AbstractArray,phi::AbstractArray)
integrate(f,x,w,j)
integrate(f::AbstractArray,x,w,j)
```

## Differential operators

In addition to the `gradient` function already discussed, the following differential operators
are defined.

```@docs
divergence(f)
symmetric_gradient(f)
ε
curl(f)
grad2curl(f)
laplacian(f)
Δ
```

A fancy way of typing the differential operators is via the nabla object defined in the library.

```@docs
(*)(::typeof(∇),f)
outer(::typeof(∇),f)
outer(f,::typeof(∇))
cross(::typeof(∇),f)
```

## Concrete Field implementations

```@docs
Homothecy
AffineMap
```

