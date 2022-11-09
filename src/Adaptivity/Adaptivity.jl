"""

"""
module Adaptivity

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

import Base: view, *
import LinearAlgebra: mul!
import Gridap.Geometry: get_grid, get_grid_topology, get_face_labeling
import Gridap.Geometry: Triangulation, is_change_possible, best_target, get_background_model
import Gridap.Geometry: move_contributions
import Gridap.CellData: change_domain

include("AdaptedDiscreteModels.jl")
include("AdaptedTriangulations.jl")
include("CompositeMeasures.jl")
include("GridTransferOperators.jl")

export AdaptivityGlue
export get_f2c_ref_cell_map, get_f2c_ref_coordinate_map

export AdaptedDiscreteModel
export get_model, get_parent, get_adaptivity_glue
export refine

export AdaptedTriangulation
export Triangulation, is_change_possible, best_target, get_adapted_model
export change_domain_c2f, change_domain, move_contributions

export CompositeMeasure, integrate

export ProjectionTransferOperator
export mul!


end # module
