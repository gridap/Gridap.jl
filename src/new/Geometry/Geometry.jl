"""

Exported names are
$(EXPORTS)
"""
module Geometry

using Test
using DocStringExtensions
using FillArrays

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs

import Gridap.ReferenceFEs: get_node_coordinates
import Gridap.ReferenceFEs: num_nodes
import Gridap.ReferenceFEs: is_affine
import Gridap.ReferenceFEs: has_straight_faces

export Triangulation
export get_cell_coordinates
export get_reffes
export get_cell_types
export num_cells
export num_cell_dims
export num_point_dims
export get_cell_reffes
export get_cell_shapefuns
export get_cell_map
export test_triangulation

export ConformingTriangulation
export get_cell_nodes
export test_conforming_triangulation

export UnstructuredGrid

include("Triangulations.jl")

include("ConformingTriangulations.jl")

include("UnstructuredGrids.jl")

end # module
