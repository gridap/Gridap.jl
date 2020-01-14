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
using Gridap.Polynomials
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
import Gridap.Arrays: reindex
import Gridap.Arrays: get_array

import Gridap.TensorValues: outer
import Gridap.TensorValues: inner
import Gridap.TensorValues: symmetic_part
import Gridap.Fields: gradient
import Gridap.Fields: grad2curl
import Base: +, - , *
import LinearAlgebra: cross
import LinearAlgebra: tr
import Base: transpose
import Base: adjoint

import Gridap.Io: to_dict
import Gridap.Io: from_dict

using Gridap.Fields: AffineMapGrad

import Gridap.ReferenceFEs: get_node_coordinates
import Gridap.ReferenceFEs: num_nodes
import Gridap.ReferenceFEs: is_affine
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

import Gridap.Fields: field_array_gradient
import Gridap.Fields: apply_lincomb
import Gridap.Fields: evaluate_field_array
import Gridap.Fields: kernel_evaluate
import Gridap.Fields: evaluate

import Gridap.Arrays: apply_kernel!

export CellField
export GenericCellField
export similar_cell_field
export test_cell_field
export convert_to_cell_field

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
export RegularityStyle
export is_oriented
export is_regular

export UnstructuredGridTopology

export Triangulation
export get_reffes
export get_cell_coordinates
export get_cell_reffes
export get_cell_shapefuns
export get_cell_map
export get_normal_vector
export test_triangulation
export restrict

export Grid
export get_cell_nodes
export test_grid
export compute_linear_grid
export compute_reference_grid

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

export DiscreteModel
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
export get_volume_triangulation
export get_face_to_cell
export get_face_to_cell_map
export test_boundary_triangulation

export GenericBoundaryTriangulation

export SkeletonPair
export jump
export mean
export SkeletonTriangulation


include("GridTopologies.jl")

include("GridTopologyMocks.jl")

include("UnstructuredGridTopologies.jl")

include("CellFields.jl")

include("Triangulations.jl")

include("Grids.jl")

include("GridMocks.jl")

include("GridPortions.jl")

include("UnstructuredGrids.jl")

include("CartesianGrids.jl")

include("FaceLabelings.jl")

include("DiscreteModels.jl")

include("DiscreteModelMocks.jl")

include("UnstructuredDiscreteModels.jl")

include("CartesianDiscreteModels.jl")

include("BoundaryTriangulations.jl")

include("GenericBoundaryTriangulations.jl")

include("SkeletonPairs.jl")

include("SkeletonTriangulations.jl")

end # module
