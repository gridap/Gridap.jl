module Fields

#using Gridap.Arrays: print_op_tree
#using Gridap.Arrays: Map
#using Gridap.Arrays: evaluate
#using Gridap.Arrays: Operation
#using Gridap.Arrays: Broadcasting
#using Gridap.Arrays: test_map
#using Gridap.Arrays: LazyArray
#using Gridap.Arrays: CachedArray
#using Gridap.Arrays: setsize!
#using Gridap.Arrays: get_array
#using Gridap.Arrays: testitem
#using Gridap.Arrays: TreeNode
#using Gridap.Arrays: similar_tree_node

using Gridap.Arrays
import Gridap.Arrays: testvalue
import Gridap.Arrays: inverse_map
import Gridap.Arrays: get_children
import Gridap.Arrays: is_zero_block
import Gridap.Arrays: testitem

using Gridap.Helpers: @abstractmethod, @notimplemented
using Gridap.Helpers: @notimplementedif, @unreachable, @check

using Gridap.Algebra: mul!
using Gridap.Algebra: fill_entries!

using Gridap.TensorValues

using LinearAlgebra: mul!, Transpose

using ForwardDiff
using FillArrays
using Test
using BlockArrays

import LinearAlgebra: det, inv, transpose, tr, cross
import LinearAlgebra: ⋅, dot

import Base: +, -, *, /
import Gridap.TensorValues: ⊗, ⊙, symmetric_part, outer

import Gridap.Arrays: IndexStyle
import Gridap.Arrays: return_cache
import Gridap.Arrays: return_type
import Gridap.Arrays: testargs
import Gridap.Arrays: return_value
import Gridap.Arrays: evaluate!
import Gridap.Arrays: lazy_map
import Gridap.Arrays: array_cache
import Gridap.Arrays: is_zero_block
# import Gridap.Arrays: uses_hash

export evaluate
export evaluate!
export return_type
export return_cache

export Field
export GenericField
export ConstantField
export FieldGradient
#export BroadcastField
export ZeroField
export MockField
export MockFieldArray
export Point
export inverse_map

export AffineMap

export gradient
export ∇
export ∇∇
export gradient_type
export push_∇
export push_∇∇

export curl
export grad2curl
export laplacian
export divergence
export Δ
export ε
export symmetric_gradient

export test_field
export test_field_array
#export test_operation_field_array
#export test_broadcast_field_array

#export mock_field

#export TransposeFieldVector
export TransposeFieldIndices
#export BroadcastOpFieldArray
#export DotOpFieldVectors
#export LinearCombinationField
#export CompositionFieldArrayField
export FieldGradientArray
#export FieldHessianArray
export BroadcastingFieldOpMap

export linear_combination
export TransposeMap
export LinearCombinationMap
export LinearCombinationField
export LinearCombinationFieldVector
export integrate
export IntegrationMap

export BlockFieldArrayCoo
export BlockFieldArrayCooMap
export similar_range

##export MatMul
#export LinCombVal
#export Integrate

include("FieldsInterfaces.jl")

include("FieldArrays.jl")

include("MockFields.jl")

include("AffineMaps.jl")

include("ApplyOptimizations.jl")

include("DiffOperators.jl")

include("AutoDiff.jl")

include("BlockFieldArrays.jl")

#include("AlgebraMaps.jl")

end
