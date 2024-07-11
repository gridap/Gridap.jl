"""
  Mesh Adaptivity for Gridap
"""
module Adaptivity

using Test
using FillArrays
using LinearAlgebra
using DataStructures

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Visualization

import Base: view, *
import Gridap.Geometry: get_grid, get_grid_topology, get_face_labeling
import Gridap.Geometry: Triangulation, is_change_possible, best_target, get_background_model
import Gridap.Geometry: move_contributions
import Gridap.CellData: change_domain

export RefinementRule
export get_cell_map, get_inverse_cell_map

export AdaptivityGlue
export get_n2o_reference_coordinate_map

export AdaptedDiscreteModel
export get_model, get_parent, get_adaptivity_glue
export is_child, is_related
export refine, coarsen, adapt

export AdaptedTriangulation
export Triangulation, is_change_possible, best_target, get_adapted_model
export change_domain, move_contributions

include("RefinementRules.jl")
include("FineToCoarseFields.jl")
include("OldToNewFields.jl")
include("FineToCoarseReferenceFEs.jl")
include("AdaptivityGlues.jl")
include("AdaptedDiscreteModels.jl")
include("AdaptedTriangulations.jl")
include("CompositeQuadratures.jl")
include("EdgeBasedRefinement.jl")
include("SimplexifyRefinement.jl")

end # module
