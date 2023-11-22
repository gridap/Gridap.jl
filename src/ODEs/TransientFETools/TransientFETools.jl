"""

The exported names are
$(EXPORTS)
"""
module TransientFETools

using Test

using DocStringExtensions

import Base: iterate
import Base: getindex
import Base: getproperty
import Base: length

using BlockArrays

using Gridap.Helpers
using Gridap.Algebra
using Gridap.CellData
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.MultiField
using Gridap.ODETools

import Gridap.Algebra: allocate_residual
import Gridap.Algebra: residual!
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: jacobian!
import Gridap.Algebra: solve

import Gridap.Fields: evaluate
import Gridap.Fields: evaluate!

import Gridap.Polynomials: get_order

import Gridap.CellData: get_data
import Gridap.CellData: get_triangulation
import Gridap.CellData: DomainStyle
import Gridap.CellData: gradient
import Gridap.CellData: ∇∇
import Gridap.CellData: change_domain

import Gridap.FESpaces: zero_free_values
import Gridap.FESpaces: get_vector_type
import Gridap.FESpaces: get_dof_value_type
import Gridap.FESpaces: has_constraints
import Gridap.FESpaces: BasisStyle
import Gridap.FESpaces: get_test
import Gridap.FESpaces: get_trial
import Gridap.FESpaces: SparseMatrixAssembler
import Gridap.FESpaces: get_algebraic_operator

using Gridap.MultiField: MultiFieldFESpace
using Gridap.MultiField: ConsecutiveMultiFieldStyle
using Gridap.MultiField: BlockMultiFieldStyle
using Gridap.MultiField: BlockSparseMatrixAssembler

import Gridap.MultiField: num_fields
import Gridap.MultiField: MultiFieldStyle
import Gridap.MultiField: ConstraintStyle

import Gridap.ODEs.ODETools: ∂t
import Gridap.ODEs.ODETools: ∂tt
import Gridap.ODEs.ODETools: ODEOperatorType
import Gridap.ODEs.ODETools: jacobians!
import Gridap.ODEs.ODETools: allocate_cache
import Gridap.ODEs.ODETools: update_cache!

include("TransientFESpaces.jl")

export TransientTrialFESpace
export TransientMultiFieldTrialFESpace
export TransientMultiFieldFESpace

export ∂t
export ∂tt

export test_transient_trial_fe_space

include("TransientCellField.jl")

export TransientCellField

include("TransientMultiFieldCellField.jl")

include("TransientFEOperators.jl")

export TransientFEOperator
# export TransientAffineFEOperator
# export TransientConstantFEOperator
# export TransientConstantMatrixFEOperator
# export TransientRungeKuttaFEOperator
# export TransientIMEXRungeKuttaFEOperator
# export TransientEXRungeKuttaFEOperator

export test_transient_fe_operator

include("ODEOperatorInterfaces.jl")

include("TransientFESolutions.jl")

export TransientFESolution

export test_transient_fe_solution
export test_transient_fe_solver

end # module TransientFETools
