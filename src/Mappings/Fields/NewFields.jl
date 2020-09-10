module NewFields

using Gridap.Helpers
using Gridap.Inference
using Gridap.Mappings
using Gridap.Arrays
using Gridap.Arrays: CachedArray
using Gridap.TensorValues: MultiValue, VectorValue, mutable, outer

using Test

using FillArrays

import Gridap.Mappings: return_type
import Gridap.Mappings: return_cache
import Gridap.Mappings: evaluate!
import Gridap.Mappings: evaluate
import Gridap.Mappings: return_caches
import Gridap.Mappings: return_types
import Gridap.Mappings: test_mapping
import Gridap.Mappings: gradient, ∇

export NewField
export ConstantField
export FunctionField
export LinearCombinationField
export MockField

# export gradient
# export ∇
# export gradients
export test_field
export Point

export evaluate_gradient!
export return_gradient_cache
export evaluate_hessian!
export return_hessian_cache

export linear_combination

export mock_field

include("FieldsInterfaces.jl")

include("MockFields.jl")

include("ConstantFields.jl")

include("FunctionFields.jl")

include("LinearCombinationFields.jl")

end
