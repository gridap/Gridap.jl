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
using Gridap.TensorValues

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

using Gridap.Fields: HomothecyGrad

import Gridap.ReferenceFEs: get_node_coordinates
import Gridap.ReferenceFEs: num_nodes
import Gridap.ReferenceFEs: is_affine
import Gridap.ReferenceFEs: has_straight_faces
import Gridap.ReferenceFEs: get_reffes
import Gridap.ReferenceFEs: get_faces
import Gridap.ReferenceFEs: get_face_vertices
import Gridap.ReferenceFEs: num_dims
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

import Gridap.Fields: field_array_gradient

export GridTopology
export compute_reffaces
export compute_cell_faces
export compute_face_vertices
export compute_isboundary_face
export test_grid_topology

export Triangulation
export get_cell_coordinates
export get_cell_type
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
export OrientationStyle
export is_oriented
export ConformityStyle
export RegularConformity
export IrregularHConformity
export IrregularPConformity
export IrregularHPConformity

export DiscreteModel
export get_node_face_owner
export get_isboundary_face
export get_face_reffe_type
export get_face_polytope_type
export get_polytopes
export get_vertex_coordinates
export get_cell_faces
export get_isboundary_node
export get_face_labeling
export get_cell_perm_indices
export get_reffes_offsets
export get_reffes_alldims
export extract_face_reffes
export replace_reffes
export test_discrete_model

export UnstructuredGrid
export UnstructuredDiscreteModel

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
export add_tag!
export add_tag_from_tags!

export CartesianDiscreteModel

include("GridTopologies.jl")

include("GridTopologyMocks.jl")

include("Triangulations.jl")

include("ConformingTriangulations.jl")

include("ConformingTrianMocks.jl")

include("UnstructuredGrids.jl")

include("CartesianGrids.jl")

include("DiscreteModels.jl")

include("FaceLabelings.jl")

include("DiscreteModelMocks.jl")

include("UnstructuredDiscreteModels.jl")

include("CartesianDiscreteModels.jl")

end # module
