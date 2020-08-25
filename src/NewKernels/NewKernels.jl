module NewKernels

using Gridap.Helpers
using Gridap.Inference

using Test

export NewKernel

export FunctionKernel
export BroadcastKernel
export CompositionKernel

export return_type
export cache
export evaluate!
export evaluate
export test_kernel
export caches
export return_types

include("Interface.jl")

include("BroadcastKernels.jl")

include("FunctionKernels.jl")

include("CompositionKernels.jl")

end
