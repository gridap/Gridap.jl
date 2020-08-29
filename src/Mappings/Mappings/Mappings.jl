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

export return_type
export return_cache
export evaluate!
export evaluate
export test_mapping
export return_caches
export return_types

export MappingArray

export apply
export apply_all

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

include("MappingArrays.jl")

end
