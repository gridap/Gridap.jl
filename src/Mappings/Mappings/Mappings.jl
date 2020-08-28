module Mappings

using Gridap.Helpers
using Gridap.Inference

using Gridap.Arrays

using Gridap.Algebra: mul!

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

# Field

include("MappingInterfaces.jl")

include("BroadcastMappings.jl")

include("FunctionMappings.jl")

include("CompositionMappings.jl")

include("AlgebraMappings.jl")

end
