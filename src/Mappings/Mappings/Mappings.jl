module Mappings

using Gridap.Helpers
using Gridap.Inference

using Gridap.Arrays

using Gridap.Algebra: mul!

using FillArrays

using Test

import Gridap.Inference: return_type

export Mapping

export FunctionMapping
export BroadcastMapping
export CompositionMapping
export MatVecMapping
export ConstantMapping

export return_type
export return_cache
export evaluate!
export evaluate
export multievaluate
export test_mapping
export return_caches
export return_types
export gradient

export MappedArray

export apply_mapping
export apply_mappings
export apply_function
export test_mapped_array
# export apply_all

export CompositionMapping

export composition

# To create the MappingArray <: AbstractArray
# using the extended Gridap API
import Gridap.Arrays: array_cache
import Gridap.Arrays: testitem
import Gridap.Arrays: uses_hash
import Gridap.Arrays: getindex!
# @santiagobadia : Probably just export in the future
import Gridap.Arrays: apply

# Field

include("MappingInterfaces.jl")

include("BroadcastMappings.jl")

include("FunctionMappings.jl")

include("CompositionMappings.jl")

include("AlgebraMappings.jl")

include("MappedArrays.jl")

include("MappingArrays.jl")

include("ConstantMappings.jl")

end
