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
- 2nd order, symmetric and traceless [`SymTracelessTensorValue`](@ref),
- 3rd order [`ThirdOrderTensorValue`](@ref),
- 4th order and symmetric [`SymFourthOrderTensorValue`](@ref).

## Generalities

The main feature of this module is that the provided types do not extend from `AbstractArray`, but from `Number`!

This allows one to work with them as if they were scalar values in broadcasted operations on arrays of `VectorValue` objects (also for `TensorValue` or `MultiValue` objects). For instance, one can perform the following manipulations:
```julia
# Assign a VectorValue to all the entries of an Array of VectorValues
A = zeros(VectorValue{2,Int}, (4,5))
v = VectorValue(12,31)
A .= v # This is possible since  VectorValue <: Number

# Broadcasting of tensor operations in arrays of TensorValues
t = TensorValue(13,41,53,17) # creates a 2x2 TensorValue
g = TensorValue(32,41,3,14) # creates another 2x2 TensorValue
B = fill(t,(1,5))
C = inner.(g,B) # inner product of g against all TensorValues in the array B
@show C
# C = [2494 2494 2494 2494 2494]
```

To create a [`::MultiValue`](@ref) tensor from components, these should be given
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
For symmetric tensor types, only the independent components should be given, see
[`SymTensorValue`](@ref), [`SymTracelessTensorValue`](@ref) and [`SymFourthOrderTensorValue`](@ref).

A `MultiValue` can be created from an `AbstractArray` of the same size. If the
`MultiValue` type has internal constraints (e.g. symmetries), ONLY the required
components are picked from the array WITHOUT CHECKING if the given array
did respect the constraints:
```julia
SymTensorValue( [1 2; 3 4] )          # -> SymTensorValue{2, Int64, 3}(1, 2, 4)
SymTensorValue( SMatrix{2}(1,2,3,4) ) # -> SymTensorValue{2, Int64, 3}(1, 3, 4)
```

`MultiValue`s can be converted to static and mutable arrays types from
`StaticArrays.jl` using `convert` and [`mutable`](@ref), respectively.

## Tensor types

The following concrete tensor types are currently implemented:

```@docs
VectorValue
TensorValue
SymTensorValue
SymTracelessTensorValue
ThirdOrderTensorValue
SymFourthOrderTensorValue
```

### Abstract tensor types

```@docs
MultiValue
AbstractSymTensorValue
```

## Interface

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
indep_components_names
change_eltype

inner
dot
double_contraction
outer
```

### Other type specific interfaces

#### For square second order tensors

```@docs
det
inv
symmetric_part
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
cross
```

#### For second order non-traceless and symmetric fourth order tensors

```@docs
one
```

