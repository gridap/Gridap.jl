"""

The exported names are
$(EXPORTS)
"""
module MultiField

using DocStringExtensions
using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields

using Gridap.Fields: ArrayBlock, MatrixBlock, VectorBlock
using Gridap.FESpaces: FEBasis, TestBasis, TrialBasis, get_cell_dof_values
using Gridap.FESpaces: SingleFieldFEBasis, TestBasis, TrialBasis
using Gridap.FESpaces: get_rows, get_cols, get_assembly_strategy
using Gridap.CellData: CellFieldAt
using Gridap.CellData: SkeletonCellFieldPair

using SparseMatricesCSR: SparseMatrixCSR

import Gridap.Fields: gradient, DIV, ∇∇

using ForwardDiff
using FillArrays
using SparseArrays
using LinearAlgebra
using BlockArrays

export num_fields
export compute_field_offsets
export restrict_to_field
export MultiFieldCellField
export MultiFieldFESpace
export MultiFieldFEFunction
export MultiFieldStyle
export ConsecutiveMultiFieldStyle
export BlockMultiFieldStyle

include("MultiFieldCellFields.jl")

include("MultiFieldFESpaces.jl")

include("MultiFieldFEFunctions.jl")

include("MultiFieldFEAutodiff.jl")

include("BlockSparseMatrixAssemblers.jl")

end # module
