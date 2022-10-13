"""

"""
module Refinement

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Visualization

import Base: view
import Gridap.Geometry: get_grid, get_grid_topology, get_face_labeling
import Gridap.Geometry: Triangulation, is_change_possible, best_target, get_glue, get_background_model
import Gridap.CellData: change_domain

include("RefinedDiscreteModels.jl")
include("RefinedTriangulations.jl")

export RefinementGlue
export get_f2c_ref_cell_map, get_f2c_ref_coordinate_map

export RefinedDiscreteModel
export get_model, get_parent, get_glue
export RefinedCartesianDiscreteModel

export RefinedTriangulation
export Triangulation, is_change_possible, best_target, get_refined_model
export change_domain_c2f, change_domain

end # module
