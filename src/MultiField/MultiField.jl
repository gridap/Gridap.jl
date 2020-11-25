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
using Gridap.Integration
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields

using Gridap.FESpaces: FEBasis, TestBasis, TrialBasis
using Gridap.Arrays: BlockArrayCooMap

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

include("MultiFieldCellFields.jl")

include("MultiFieldFESpaces.jl")

include("MultiFieldFEFunctions.jl")

end # module
