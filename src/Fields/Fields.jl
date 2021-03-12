module Fields

using Gridap.Arrays
import Gridap.Arrays: testvalue
import Gridap.Arrays: inverse_map
import Gridap.Arrays: get_children
import Gridap.Arrays: is_zero_block
import Gridap.Arrays: testitem
using Gridap.Arrays: BlockArrayCooMap

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
using StaticArrays

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

export evaluate
export evaluate!
export return_type
export return_cache

export Field
export GenericField
export ConstantField
export constant_field
export FieldGradient
export FieldGradientArray
export ZeroField
export MockField
export MockFieldArray
export Point
export inverse_map

export AffineMap
export affine_map

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

export linear_combination
export integrate
export IntegrationMap

include("FieldsInterfaces.jl")

include("FieldArrays.jl")

include("MockFields.jl")

include("AffineMaps.jl")

include("ApplyOptimizations.jl")

include("DiffOperators.jl")

include("AutoDiff.jl")

include("BlockFieldArrays.jl")

end
