"""

"""
module Refinement

using FillArrays
using LinearAlgebra
using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Visualization

import Base: view
import LinearAlgebra: mul!
import Gridap.Geometry: get_grid, get_grid_topology, get_face_labeling
import Gridap.Geometry: Triangulation, is_change_possible, best_target, get_background_model
import Gridap.CellData: change_domain

include("RefinedDiscreteModels.jl")
include("RefinedTriangulations.jl")
include("GridTransferOperators.jl")

export RefinementGlue
export get_f2c_ref_cell_map, get_f2c_ref_coordinate_map

export RefinedDiscreteModel
export get_model, get_parent, get_refinement_glue
export refine

export RefinedTriangulation
export Triangulation, is_change_possible, best_target, get_refined_model
export change_domain_c2f, change_domain

export ProjectionTransferOperator, RefinementTransferMap
export mul!


end # module
