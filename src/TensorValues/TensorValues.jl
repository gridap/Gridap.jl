"""
Immutable tensor types for Gridap.
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
