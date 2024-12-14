"""
Immutable tensor types for Gridap. The currently implemented tensor types are
- 1st order [`VectorValue`](@ref),
- 2nd order [`TensorValue`](@ref),
- 2nd order and symmetric [`SymTensorValue`](@ref),
- 2nd order, symmetric and traceless [`SymTracelessTensorValue`](@ref),
- 3rd order [`ThirdOrderTensorValue`](@ref),
- 4th order and symmetric [`SymFourthOrderTensorValue`](@ref).

Example usage:
```julia
# create a 2D vector from components
v = VectorValue(12,31)

# Assign a VectorValue to all the entries of an Array of VectorValues
A = zeros(VectorValue{2,Int}, (4,5))
A .= v # This is possible since  VectorValue <: Number

using StaticArrays
# create 2x2 tensor from component tuple
t = TensorValue( (1, 2, 3, 4) )
# conversion to StaticArrays type
ts= convert(SMatrix{2,2,Int}, t)
@show ts
# 2×2 SMatrix{2, 2, Int64, 4} with indices SOneTo(2)×SOneTo(2):
#  1  3
#  2  4
t2[1,2] == t[1,2] == 3 # true

# conversion from Array or StaticArray types, symmetric tensor types only store required components
SymTensorValue( [1 2; 3 4] )          # SymTensorValue{2, Int64, 3}(1, 2, 4)
SymTensorValue( SMatrix{2}(1,2,3,4) ) # SymTensorValue{2, Int64, 3}(1, 3, 4)
```

See the official documentation for more details.
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

export MultiValue
export VectorValue
export TensorValue
export AbstractSymTensorValue
export SymTensorValue
export SymTracelessTensorValue
export QTensorValue
export SymFourthOrderTensorValue
export ThirdOrderTensorValue

export inner, outer, meas
export mutable
export Mutable
export symmetric_part
export n_components
export num_components
export num_indep_components
export change_eltype
export diagonal_tensor
export ⊙
export ⊗
export ⋅¹
export ⋅²
export double_contraction
export data_index
export indep_comp_getindex
export indep_components_names

import Base: show
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

include("SymFourthOrderTensorValueTypes.jl")

include("ThirdOrderTensorValueTypes.jl")

include("Indexing.jl")

include("Operations.jl")

include("Reinterpret.jl")

end # module
