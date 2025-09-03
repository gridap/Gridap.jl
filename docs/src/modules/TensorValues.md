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

## Interface and operations

The tensor types implement methods for the following `Base` functions: `getindex`, `length`, `size`, `rand`, `zero`, `real`, `imag` and `conj`.

`one` is also implemented in particular cases: it is defined for second
and fourth order tensors. For second order, it returns the identity tensor `δij`,
except `SymTracelessTensorValue` that does not implement `one`. For fourth order symmetric tensors, see [`one`](@ref).

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
