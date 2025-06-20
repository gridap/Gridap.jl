"""
This module provide tensor types of supertype [`MultiValue`](@ref) that have
very similar design and API to `S/MArray` from StaticArrays.jl, but with
efficient storage for tensors with dependent components (e.g. symmetry). The
main feature of this module is that `MultiValue` doesn't subtype `AbstractArray`,
but `Number`!

This allows one to work with them as if they were scalar values in broadcasted
operations on arrays of `VectorValue` objects (also for `TensorValue` or any
`<:MultiValue` objects). For instance, one can perform the following manipulations:
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
SymTensorValue( [1 3; 2 4] )          # -> SymTensorValue{2, Int64, 3}(1, 3, 4)
SymTensorValue( SMatrix{2}(1,2,3,4) ) # -> SymTensorValue{2, Int64, 3}(1, 3, 4)
```

`MultiValue`s can be converted to static and mutable arrays types from
`StaticArrays.jl` using `convert` and [`mutable`](@ref), respectively.

$(public_names_in_md(@__MODULE__; change_link=Dict(
  :QTensorValue  => "SymTracelessTensorValue",
  :×  => "cross",
  :⊗  => "outer",
  :⊙  => "inner",
  :⋅  => "dot",
  :⋅¹ => "dot",
  :⋅² => "double_contraction",
)))
"""
module TensorValues

using DocStringExtensions

using StaticArrays
using StaticArrays: SVector, MVector, SMatrix, MMatrix, SArray, MArray
using Base: @propagate_inbounds, @pure
using Gridap.Helpers
using Gridap.Arrays
using LinearAlgebra
using Random
using ForwardDiff

export MultiValue
export VectorValue
export TensorValue
export AbstractSymTensorValue
export SymTensorValue
export SymTracelessTensorValue
export SkewSymTensorValue
export QTensorValue
export SymFourthOrderTensorValue
export ThirdOrderTensorValue

export inner, outer, meas
export mutable
export Mutable
export symmetric_part
export skew_symmetric_part
export num_components
export num_indep_components
export change_eltype
export diagonal_tensor
export ⊙
export ⊗
export ⋅¹
export ⋅²
export double_contraction
export congruent_prod
export indep_comp_getindex
export indep_components_names
export component_basis
export representatives_of_componentbasis_dual

import Base: show
import Base: promote_rule
import Base: zero, one
import Base: +, -, *, /, \, ==, ≈, isless
import Base: conj, real, imag
import Base: sum, maximum, minimum
import Base: getindex, iterate, eachindex, lastindex
import Base: size, length, eltype
import Base: reinterpret
import Base: convert
import Base: CartesianIndices
import Base: LinearIndices
import Base: adjoint
import Base: transpose
import Base: rand

import LinearAlgebra: det, inv, tr, cross, dot, norm
# Reexport from LinearAlgebra (just for convenience)
export det, inv, tr, cross, dot, norm, ×, ⋅

import Gridap.Arrays: get_array

include("MultiValueTypes.jl")

include("VectorValueTypes.jl")

include("TensorValueTypes.jl")

include("SymTensorValueTypes.jl")

include("SymTracelessTensorValueTypes.jl")

include("SkewSymTensorValueTypes.jl")

include("SymFourthOrderTensorValueTypes.jl")

include("ThirdOrderTensorValueTypes.jl")

include("Indexing.jl")

include("Operations.jl")

include("Reinterpret.jl")

end # module
