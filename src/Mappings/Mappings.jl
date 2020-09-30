module Mappings

using Gridap.Helpers
using Gridap.Inference
using Gridap.Arrays
using Gridap.Algebra: mul!
using FillArrays
using Test

import Gridap.Inference: return_type

# Mapping interface

export Mapping
export return_type
export return_cache
export evaluate!
export evaluate
export test_mapping
export return_caches
export return_types
import Gridap.Arrays: testitem
import Gridap.Arrays: getindex!
import Gridap.Arrays: array_cache
import Gridap.Arrays: testitem
import Gridap.Arrays: uses_hash
import Gridap.Arrays: getindex!
import Gridap.Arrays: IndexStyle
# import Gridap.Arrays: apply

# MappedArray

export MappedArray
export apply
export test_mapped_array

# BroadcastMapping

export BroadcastMapping

# OperationMappings

export OperationMapping
export operation
export composition
export Operation

# export derivative
# export gradient
# export ∇

# export MatVecMapping

# # To create the MappingArray <: AbstractArray
# # using the extended Gridap API
# import Gridap.Arrays: array_cache
# import Gridap.Arrays: testitem
# import Gridap.Arrays: uses_hash
# import Gridap.Arrays: getindex!
# @santiagobadia : Probably just export in the future
# import Gridap.Arrays: apply

# Field

using Gridap.TensorValues

using ForwardDiff

import LinearAlgebra: det, inv, transpose
import Base: +, -, *, /
import LinearAlgebra: ⋅

import Base: promote_type
using LinearAlgebra: mul!, Transpose

export Field
export GenericField
export Gradient
export Hessian
export ConstantField
export FunctionField
export BroadcastField
export ZeroField
export MockField
export Point

export TransposeFieldVector
export BroadcastOpFieldArray
export CompositionFieldArrayField
export DotOpFieldVectors

export evaluate_gradient!
export return_gradient_type
export return_gradient_cache
export evaluate_hessian!
export return_hessian_cache

export gradient
export ∇
export hessian

export test_field
export test_field_array
export test_operation_field_array
export test_broadcast_field_array

export mock_field

export linear_combination
# export MockBasis
# export LinearCombinationField
# export OtherMockField
# export OtherMockBasis

# export GenericFieldArray

# export FieldHessian
# export FieldGradientArray
# export FieldHessianArray


# # export gradients



# export field_composition
# export field_array_composition

# export mock_field

include("MappingInterfaces.jl")

include("MappedArrays.jl")

include("FieldsInterfaces.jl")

include("FieldArrays.jl")

include("MockFields.jl")

include("AutoDiff.jl")

# include("FunctionMappings.jl")


# include("AlgebraMappings.jl")


# include("MappingArrays.jl")

# include("ConstantMappings.jl")

# @santiagobadia : To be decided what to do here
# include("MappingGradients.jl")

end
