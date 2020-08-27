module NewFields

using Gridap.Helpers
using Gridap.Inference
using Gridap.NewKernels
using Gridap.Arrays
using Gridap.TensorValues: VectorValue, mutable, outer

using Gridap.Arrays: CachedArray

using Test

using FillArrays

import Gridap.NewKernels: return_type
import Gridap.NewKernels: return_cache
import Gridap.NewKernels: evaluate!
import Gridap.NewKernels: evaluate
import Gridap.NewKernels: return_caches
import Gridap.NewKernels: return_types
import Gridap.NewKernels: test_kernel

export NewField
export ConstantField
export FunctionField
export LinearCombinationField

export gradient
export ∇
export gradients
export test_field
export Point
export return_gradient_type

export evaluate_gradient!
export return_gradient_cache
export evaluate_hessian!
export return_hessian_cache

export linear_combination

include("Interfaces.jl")

include("MockFields.jl")

include("ConstantFields.jl")

include("FunctionFields.jl")

include("LinearCombinationFields.jl")

end
