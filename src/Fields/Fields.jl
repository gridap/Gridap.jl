module Fields

using Gridap.Arrays
import Gridap.Arrays: testvalue
import Gridap.Arrays: inverse_map
import Gridap.Arrays: get_children
import Gridap.Arrays: testitem

using Gridap.Helpers: @abstractmethod, @notimplemented
using Gridap.Helpers: @notimplementedif, @unreachable, @check
using Gridap.Helpers: tfill

using Gridap.Algebra: mul!

using Gridap.TensorValues
using Gridap.Algebra

using LinearAlgebra: mul!, Transpose, diag

using ForwardDiff
using FillArrays
using NLsolve
using Test
using StaticArrays
using LinearAlgebra

import LinearAlgebra: det, inv, transpose, tr, cross
import LinearAlgebra: ⋅, dot

import Base: +, -, *, /
import Gridap.TensorValues: ⊗, ⊙, symmetric_part, outer, meas

import Gridap.Arrays: IndexStyle
import Gridap.Arrays: return_cache
import Gridap.Arrays: return_type
import Gridap.Arrays: testargs
import Gridap.Arrays: return_value
import Gridap.Arrays: evaluate!
import Gridap.Arrays: lazy_map
import Gridap.Arrays: array_cache

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
export pinvJt

export curl
export grad2curl
export laplacian
export divergence
export DIV
export Δ
export ε
export symmetric_gradient

export test_field
export test_field_array

export linear_combination
export integrate
export IntegrationMap

export ArrayBlock
export VectorBlock
export MatrixBlock
export BlockMap
export ArrayBlockView
export VectorBlockView
export MatrixBlockView

export VoidField
export VoidFieldMap
export VoidBasis
export VoidBasisMap

export DensifyInnerMostBlockLevelMap

include("FieldsInterfaces.jl")

include("FieldArrays.jl")

include("MockFields.jl")

include("AffineMaps.jl")

include("ApplyOptimizations.jl")

include("DiffOperators.jl")

include("AutoDiff.jl")

include("ArrayBlocks.jl")

include("InverseFields.jl")

include("DensifyInnerMostBlockLevelMaps.jl")

end
