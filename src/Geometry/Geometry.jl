"""

Exported names are
$(EXPORTS)
"""
module Geometry

using Test
using DocStringExtensions
using FillArrays

using LinearAlgebra: ⋅

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.TensorValues
using Gridap.Io

using Gridap.ReferenceFEs: _num_faces
using Gridap.ReferenceFEs: _num_facets
using Gridap.ReferenceFEs: _num_vertices
using Gridap.ReferenceFEs: _num_edges
using Gridap.ReferenceFEs: _get_facedims
using Gridap.ReferenceFEs: _get_offsets
using Gridap.ReferenceFEs: _get_offset
using Gridap.ReferenceFEs: _find_unique_with_indices

import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Arrays: get_array
import Gridap.Arrays: lazy_append
import Gridap.Arrays: lazy_map

import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!
import Gridap.Arrays: get_children

import Gridap.Io: to_dict
import Gridap.Io: from_dict

import Gridap.ReferenceFEs: ReferenceFE
import Gridap.ReferenceFEs: get_node_coordinates
import Gridap.ReferenceFEs: num_nodes
import Gridap.ReferenceFEs: is_first_order
import Gridap.ReferenceFEs: get_faces
import Gridap.ReferenceFEs: get_face_vertices
import Gridap.ReferenceFEs: get_dimranges
import Gridap.ReferenceFEs: get_dimrange
import Gridap.ReferenceFEs: num_faces
import Gridap.ReferenceFEs: num_facets
import Gridap.ReferenceFEs: num_edges
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: get_facedims
import Gridap.ReferenceFEs: get_offsets
import Gridap.ReferenceFEs: get_offset
import Gridap.ReferenceFEs: get_vertex_coordinates
import Gridap.ReferenceFEs: get_face_nodes
import Gridap.ReferenceFEs: get_face_own_nodes
import Gridap.ReferenceFEs: get_polytope
import Gridap.ReferenceFEs: get_vertex_node
import Gridap.ReferenceFEs: is_simplex
import Gridap.ReferenceFEs: is_n_cube
import Gridap.ReferenceFEs: get_reffaces
import Gridap.ReferenceFEs: get_face_type
import Gridap.ReferenceFEs: num_dims
import Gridap.ReferenceFEs: num_cell_dims
import Gridap.ReferenceFEs: num_point_dims
import Gridap.ReferenceFEs: simplexify
import Gridap.ReferenceFEs: get_facet_normal

export GridTopology
export num_cells
export get_polytopes
export get_cell_type
export get_cell_vertices
export get_reffaces_offsets
export compute_reffaces
export compute_cell_faces
export compute_face_vertices
export compute_isboundary_face
export get_cell_permutations
export compute_cell_permutations
export test_grid_topology
export get_cell_faces
export get_isboundary_face
export OrientationStyle
export Oriented
export NonOriented
export RegularityStyle
export Regular
export Irregular
export is_oriented
export is_regular
export expand_cell_data
export compress_cell_data

export UnstructuredGridTopology

export Triangulation
export TriangulationStyle
export BackgroundTriangulation
export SubTriangulation
export get_reffes
export get_cell_coordinates
export get_cell_ref_coordinates
export get_cell_reffe
export get_cell_shapefuns
export get_facet_normal
export test_triangulation
export get_cell_id
export get_cell_map
export get_background_triangulation
export get_cell_ref_map
export have_compatible_domains

export Grid
export get_cell_nodes
export test_grid
export compute_linear_grid
export compute_reference_grid

export TriangulationPortion
export RestrictedTriangulation
export GridPortion
export UnstructuredGrid

export CartesianGrid
export CartesianDescriptor
export get_cartesian_descriptor

export FaceLabeling
export num_tags
export num_entities
export get_face_entity
export get_tag_entities
export get_tag_name
export get_tag_from_name
export get_tags_from_names
export add_tag!
export add_tag_from_tags!
export get_face_mask
export get_face_tag
export get_face_tag_index
export get_cell_entity

export DiscreteModel
export DiscreteModelFromFile
export get_grid
export get_grid_topology
export get_face_labeling
export test_discrete_model
export compute_face_nodes
export compute_face_own_nodes
export compute_vertex_node
export get_node_face_owner
export compute_node_face_owner
export get_triangulation

export UnstructuredDiscreteModel
export CartesianDiscreteModel

export BoundaryTriangulation
#export get_volume_triangulation
#export get_face_to_cell
#export get_face_to_lface
#export get_face_to_cell_map
#export get_face_to_face
#export get_cell_around
#export test_boundary_triangulation

#export GenericBoundaryTriangulation

export DiscreteModelPortion

export SkeletonPair
export SkeletonTriangulation
export InterfaceTriangulation
#export get_left_boundary
#export get_right_boundary

export RestrictedDiscreteModel
export get_parent_model

export AppendedTriangulation

export GridMock

include("Triangulations.jl")

include("Grids.jl")

include("GridMocks.jl")

include("UnstructuredGrids.jl")

include("CartesianGrids.jl")

include("GridTopologies.jl")

include("GridTopologyMocks.jl")

include("UnstructuredGridTopologies.jl")

include("FaceLabelings.jl")

include("DiscreteModels.jl")

include("DiscreteModelMocks.jl")

include("UnstructuredDiscreteModels.jl")

include("CartesianDiscreteModels.jl")

include("RestrictedTriangulations.jl")

include("BoundaryTriangulations.jl")

include("SkeletonTriangulations.jl")

include("GridPortions.jl")

include("DiscreteModelPortions.jl")

include("RestrictedDiscreteModels.jl")

include("AppendedTriangulations.jl")

end # module
