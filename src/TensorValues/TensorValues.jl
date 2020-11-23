"""
This module provides concrete implementations of `Number` that represent
1st, 2nd and general order tensors.

## Why

The main feature of this module is that the provided types do not extend from `AbstractArray`, but from `Number`!

This allows one to work with them as if they were scalar values in broadcasted operations on arrays of `VectorValue` objects (also for `TensorValue` or `MultiValue` objects). For instance, one can perform the following manipulations:
```julia
# Assing a VectorValue to all the entries of an Array of VectorValues
A = zeros(VectorValue{2,Int}, (4,5))
v = VectorValue(12,31)
A .= v # This is posible since  VectorValue <: Number

# Broadcasting of tensor operations in arrays of TensorValues
t = TensorValue(13,41,53,17) # creates a 2x2 TensorValue
g = TensorValue(32,41,3,14) # creates another 2x2 TensorValue
B = fill(t,(1,5))
C = inner.(g,B) # inner product of g against all TensorValues in the array B
@show C
# C = [2494 2494 2494 2494 2494]
```

The exported names are:

$(EXPORTS)
"""
module TensorValues

using DocStringExtensions

using StaticArrays
using StaticArrays: SVector, MVector, SMatrix, MMatrix, SArray, MArray
using Base: @propagate_inbounds, @pure
using Gridap.Helpers
using Gridap.Arrays
using LinearAlgebra

export MultiValue
export VectorValue
export TensorValue
export SymTensorValue
export SymFourthOrderTensorValue
export ThirdOrderTensorValue

export inner, outer, meas
export mutable
export Mutable
export symmetric_part
export n_components
export num_components
export change_eltype
export diagonal_tensor
export ⊙
export ⊗
export ⋅¹
export ⋅²
export double_contraction

import Base: show
import Base: zero, one
import Base: +, -, *, /, \, ==, ≈, isless
import Base: conj
import Base: sum, maximum, minimum
import Base: getindex, iterate, eachindex
import Base: size, length, eltype
import Base: reinterpret
import Base: convert
import Base: CartesianIndices
import Base: LinearIndices
import Base: adjoint
import Base: transpose
#import Base: :

import LinearAlgebra: det, inv, tr, cross, dot, norm
# Reexport from LinearAlgebra (just for convenience)
export det, inv, tr, cross, dot, norm, ×, ⋅

import Gridap.Arrays: get_array

include("MultiValueTypes.jl")

include("VectorValueTypes.jl")

include("TensorValueTypes.jl")

include("SymTensorValueTypes.jl")

include("SymFourthOrderTensorValueTypes.jl")

include("ThirdOrderTensorValueTypes.jl")

include("Indexing.jl")

include("Operations.jl")

include("Reinterpret.jl")

end # module
