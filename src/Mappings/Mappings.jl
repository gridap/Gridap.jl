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
# import Gridap.Arrays: apply

# MappedArray

export MappedArray
export apply
export apply
export test_mapped_array

# BroadcastMapping

export BroadcastMapping

# OperationMappings

export OperationMapping
export operation

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

include("MappingInterfaces.jl")

include("MappedArrays.jl")

include("BroadcastMappings.jl")

include("OperationMappings.jl")

# include("FunctionMappings.jl")


# include("AlgebraMappings.jl")


# include("MappingArrays.jl")

# include("ConstantMappings.jl")

# @santiagobadia : To be decided what to do here
# include("MappingGradients.jl")

end
