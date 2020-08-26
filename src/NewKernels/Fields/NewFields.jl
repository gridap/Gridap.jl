module NewFields

using Gridap.Helpers
using Gridap.Inference
using Gridap.NewKernels
using Gridap.Arrays
using Gridap.TensorValues: VectorValue, mutable, outer

using Gridap.Arrays: CachedArray

using Test

import Gridap.NewKernels: return_type
import Gridap.NewKernels: return_cache
import Gridap.NewKernels: evaluate!
import Gridap.NewKernels: evaluate
import Gridap.NewKernels: return_caches
import Gridap.NewKernels: return_types
import Gridap.NewKernels: test_kernel

export NewField

export gradient
export gradients
export test_field
export Point
export return_gradient_type

export evaluate_gradient!
export return_gradient_cache
export evaluate_hessian!
export return_hessian_cache

include("Interfaces.jl")

include("MockFields.jl")

end
