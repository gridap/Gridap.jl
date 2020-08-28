module NewFields

using Gridap.Helpers
using Gridap.Inference
using Gridap.Mappings
using Gridap.Arrays
using Gridap.TensorValues: VectorValue, mutable, outer

using Gridap.Arrays: CachedArray

using Test

using FillArrays

import Gridap.Mappings: return_type
import Gridap.Mappings: return_cache
import Gridap.Mappings: evaluate!
import Gridap.Mappings: evaluate
import Gridap.Mappings: return_caches
import Gridap.Mappings: return_types
import Gridap.Mappings: test_mapping



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
