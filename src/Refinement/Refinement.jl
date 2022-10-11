"""

"""
module Refinement

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.Visualization

import Gridap.Geometry: get_glue

include("RefinedDiscreteModels.jl")

export RefinementGlue
export get_f2c_ref_cell_map, get_f2c_ref_coordinate_map

export RefinedDiscreteModel
export get_model, get_parent, get_glue
export change_domain_c2f

end # module
