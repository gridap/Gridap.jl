module Inverse

using Gridap, ChainRulesCore, ForwardDiff
import LinearAlgebra

import Gridap.FESpaces: FEFunction, FESpace, TrialFESpace
import Gridap.FESpaces: assemble_vector, assemble_matrix
import Gridap.FESpaces: num_free_dofs, num_dirichlet_dofs
import Gridap.FESpaces: get_cell_dof_ids, get_cell_dof_values
import Gridap.Arrays: lazy_map, get_array, ∑
import Gridap.Algebra: solve!
import Gridap.TensorValues: VectorValue
import Gridap.Fields: Point
import Gridap.CellData: ∫, ∇, ⋅, ⊙, ×
import Gridap.Algebra: numerical_setup, numerical_setup!, symbolic_setup
import Gridap.Algebra: SymbolicSetup, NumericalSetup, LinearSolver, LUSolver
import Gridap.Arrays: return_cache, array_cache, getindex!, evaluate, evaluate!
import Gridap.CellData: _point_to_cell!, get_cell_points
import Gridap.CellData: make_inverse_table, _point_to_cell_cache, get_data
import Gridap.CellData: KDTreeSearch, Broadcasting, DomainContribution, Reindex
import Gridap.CellData: GenericCellField, CellField, Measure
import Gridap.FESpaces: SingleFieldFESpace, ZeroMeanFESpace, DirichletFESpace
import Gridap.FESpaces: FESpaceWithLinearConstraints, HomogeneousTrialFESpace
import Gridap.FESpaces: SparseMatrixAssembler, assemble_matrix_and_vector
import Gridap.FESpaces: assemble_vector!, collect_cell_matrix_and_vector
import Gridap.FESpaces: zero_free_values, get_free_dof_values, get_dirichlet_dof_values
import Gridap.FESpaces: get_fe_basis, get_trial_fe_basis
import Gridap.FESpaces: get_fe_dof_basis, get_vector_type, _DOF_to_dof
import Gridap.Fields: inverse_map
import Gridap.Geometry: get_cell_map
import Gridap.Helpers: @notimplemented, @abstractmethod, @unreachable
import Gridap.ReferenceFEs: LagrangianDofBasis, MomentBasedDofBasis
import Gridap.MultiField: MultiFieldFESpace, MultiFieldStyle
import Gridap.MultiField: StridedMultiFieldStyle, ConsecutiveMultiFieldStyle
import Gridap.MultiField: num_fields

import ChainRulesCore: Tangent, NoTangent, unthunk

import FillArrays: Fill

import SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, sparse, spzeros, droptol!

import LinearAlgebra: Diagonal, mul!, norm1, norm2, norm_sqr

import Random: rand!

export norm1, norm2, norm_sqr

include("FEObservationOperators.jl")
export FEObservationOperator
export MultiFieldDataMisfitCalculator
export filter_observation_values

end
