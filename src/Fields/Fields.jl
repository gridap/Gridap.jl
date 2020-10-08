module Fields

using Gridap.Arrays: Mapping
using Gridap.Arrays: evaluate
using Gridap.Arrays: Operation
using Gridap.Arrays: Broadcasting
using Gridap.Arrays: test_mapping
using Gridap.Arrays: LazyArray
using Gridap.Arrays: test_lazy_array
using Gridap.Arrays: CachedArray
using Gridap.Arrays: setsize!
using Gridap.Arrays: get_array

using Gridap.Helpers: @abstractmethod, @notimplemented
using Gridap.Helpers: @notimplementedif, @unreachable

using Gridap.Algebra: mul!

using Gridap.TensorValues

using LinearAlgebra: mul!, Transpose

using ForwardDiff
using FillArrays
using Test

import LinearAlgebra: det, inv, transpose
import LinearAlgebra: ⋅

import Base: +, -, *, /

import Gridap.Inference: return_type

import Gridap.Arrays: IndexStyle
import Gridap.Arrays: return_cache
import Gridap.Arrays: return_type
import Gridap.Arrays: evaluate!
import Gridap.Arrays: lazy_map
import Gridap.Arrays: array_cache
# import Gridap.Arrays: uses_hash

export evaluate!
export return_type
export return_cache

export Field
export GenericField
export FieldGradient
export FieldHessian
export BroadcastField
export ZeroField
export MockField
export Point

export evaluate_gradient!
export return_gradient_type
export return_gradient_cache
export evaluate_hessian!
export return_hessian_cache
export return_hessian_type

export gradient
export ∇
export hessian

export test_field
export test_field_array
export test_operation_field_array
export test_broadcast_field_array

export mock_field

export TransposeFieldVector
export TransposeFieldIndices
export BroadcastOpFieldArray
export DotOpFieldVectors
export LinearCombinationField
export CompositionFieldArrayField
export FieldGradientArray
export FieldHessianArray

export linear_combination
export integrate

export MatMul
export LinCombVal
export Integrate

include("FieldsInterfaces.jl")

include("MockFields.jl")

include("FieldArrays.jl")

include("ApplyOptimizations.jl")

include("AutoDiff.jl")

include("AlgebraMappings.jl")

end
