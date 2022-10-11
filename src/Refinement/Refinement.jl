"""

"""
module Refinement

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Visualization

import Gridap.Geometry: get_glue
import Gridap.Geometry: Triangulation, is_change_possible, best_target

include("RefinedDiscreteModels.jl")
include("RefinedTriangulations.jl")

export RefinementGlue
export get_f2c_ref_cell_map, get_f2c_ref_coordinate_map

export RefinedDiscreteModel
export get_model, get_parent, get_glue
export change_domain_c2f

export RefinedTriangulation
export Triangulation, is_change_possible, best_target

end # module
