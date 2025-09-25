# Gridap.TensorValues

```@meta
CurrentModule = Gridap.TensorValues
```

This module provides the abstract interface `MultiValue` representing tensors
that are also `Number`s, along with concrete implementations for the following
tensors:
- 1st order [`VectorValue`](@ref),
- 2nd order [`TensorValue`](@ref),
- 2nd order and symmetric [`SymTensorValue`](@ref),
- 2nd order and skew symmetric [`SkewSymTensorValue`](@ref),
- 2nd order, symmetric and traceless [`SymTracelessTensorValue`](@ref),
- 3rd order [`ThirdOrderTensorValue`](@ref),
- 4th order and symmetric [`SymFourthOrderTensorValue`](@ref).

## Summary

```@docs
TensorValues
```

## Construction and conversion

To create a [`MultiValue`](@ref) tensor from components, these should be given
as separate arguments or all gathered in a `tuple`. The order of the arguments
is the order of the linearized Cartesian indices of the corresponding array
(order of the `Base.LinearIndices` indices):
```julia
using StaticArrays
t = TensorValue( (1, 2, 3, 4) )
ts= convert(SMatrix{2,2,Int}, t)
@show ts
# 2×2 SMatrix{2, 2, Int64, 4} with indices SOneTo(2)×SOneTo(2):
#  1  3
#  2  4
t2[1,2] == t[1,2] == 3 # true
```
For tensor types with symmetry, only the independent components should be given,
see [`SymTensorValue`](@ref), [`SkewSymTensorValue`](@ref),
[`SymTracelessTensorValue`](@ref) and [`SymFourthOrderTensorValue`](@ref).

A `MultiValue` can be created from an `AbstractArray` of the same size. If the
`MultiValue` type has internal constraints (e.g. symmetries), ONLY the required
components are picked from the array WITHOUT CHECKING if the given array
did respect the constraints:
```julia
s1 = SymTensorValue( [1 2; 3 4] )          # -> SymTensorValue{2, Int64, 3}(1, 2, 4)
s2 = SymTensorValue( SMatrix{2}(1,2,3,4) ) # -> SymTensorValue{2, Int64, 3}(1, 3, 4)
s1 != s2 # true
s3 = SymTensorValue( (1,3,4) )             # -> SymTensorValue{2, Int64, 3}(1, 3, 4)
s4 = SymTensorValue(1, 3, 4)               # -> SymTensorValue{2, Int64, 3}(1, 3, 4)
s2 === s3 === s4 # true
```

`MultiValue`s can be converted to static and mutable arrays types from
`StaticArrays.jl` using constructors and `convert`, and also [`mutable`](@ref) for
`MArray`s only. They can also be converted to julia `Array` using the
constructor (similarly to StaticArrays).


```julia
v = VectorValue(1,2,3)
sv1 = SVector(v)
sv2 = SVector{3,Float64}(v)
sv3 = convert(SVector{3,Float64}, v)
sv1 == sv2 == sv3 # true
sv2 === sv3 # true

mv1 = MVector(v)
mv2 = convert(MVector{3,Float64}, v)
mv3 = mutable(v)
mv1 == mv2 == mv3 # true

Vector(v)          # [1, 2, 3]
Vector{Float64}(v) # [1.0, 2.0, 3.0]
Array(v)           # [1, 2, 3]
```

!!! warning
    These conversions must loose the information on internal constraints of the
    components (symmetries, etc.) because `Base.Array` and `StaticArray`s do not
    have abstractions for this.

## Tensor types

The following concrete tensor types are currently implemented:

```@docs
VectorValue
TensorValue
SymTensorValue
SkewSymTensorValue
SymTracelessTensorValue
ThirdOrderTensorValue
SymFourthOrderTensorValue
```

### Abstract tensor types

```@docs
MultiValue
AbstractSymTensorValue
```

## Indexing

```@docs
getindex(::MultiValue, ::Integer)
```

## Interface and operations

The tensor types implement methods for the following `Base` functions:
`length`, `size`, `rand`, `zero`, `real`, `imag` and `conj`.

`one` is also implemented in particular cases: it is defined for second
and fourth order tensors. For second order, it returns the identity tensor `δij`,
except `SymTracelessTensorValue` that does not implement `one`. For fourth
order symmetric tensors, see [`one`](@ref).

Additionally, the tensor types expose the following interface:

```@docs
num_components
mutable
Mutable
num_indep_components
indep_comp_getindex
component_basis
representatives_of_componentbasis_dual
indep_components_names
change_eltype

inner
dot(::MultiValue,::MultiValue)
double_contraction
outer(::MultiValue,::MultiValue)
```

### Other type specific interfaces

#### For square second order tensors

```@docs
det
inv
eigen
symmetric_part
skew_symmetric_part
congruent_prod
```

#### For first order tensors

```@docs
diagonal_tensor
```

#### For second and third order tensors

```@docs
tr
```

#### For first and second order tensors

```@docs
norm
meas
```

#### For `VectorValue` of length 2 and 3

```@docs
cross(::VectorValue,::VectorValue)
```

#### For second order non-traceless and symmetric fourth order tensors

```@docs
one
```

##### Deprecated

```@docs
data_index
n_components
```
