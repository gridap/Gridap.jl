module NewKernels

using Gridap.Helpers
using Gridap.Inference

using Gridap.Arrays: CachedArray

using Test

import Gridap.Inference: return_type

export NewKernel

export FunctionKernel
export BroadcastKernel
export CompositionKernel

export return_type
export return_cache
export evaluate!
export evaluate
export test_kernel
export return_caches
export return_types

# Field

include("Interfaces.jl")

include("BroadcastKernels.jl")

include("FunctionKernels.jl")

include("CompositionKernels.jl")

end
