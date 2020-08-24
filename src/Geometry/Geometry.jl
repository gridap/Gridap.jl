"""
The Geometry module contains data structures for space discretization.

## Mesh based finite element analysis

Mesh based finite element analysis requires multiple pieces of information:
* Information about the face structure. For instance which `vertices` belong to a given edge.
* A tagging system that allows to assign tags to parts of the space. For instance
to assign Dirichlet conditions on some parts of the boundary and von Neumann conditions on other parts.
to different boundary conditions.
* Information about the isomorphism between reference space and the physical space. This isomorphism
is encoded via the `nodes`. For linear meshes `nodes` coincide with `vertices`. The latter are by definition the "corners" of the cells.

In `Gridap` these informations are organized into the following types:

* A [`Triangulation`](@ref) contains `cells` and `nodes`. It does not know about `vertices` or the face structure.
* A [`GridTopology`](@ref) contains information about the `vertices` and the face structure.
* A [`FaceLabeling`](@ref) records which faces carry which tags.
* A [`DiscreteModel`](@ref) contains all of the above: [`Triangulation`](@ref), [`GridTopology`](@ref) and [`FaceLabeling`](@ref)

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
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.TensorValues
using Gridap.Io
using Gridap.Integration
using Gridap.CellData

using Gridap.ReferenceFEs: _num_faces
using Gridap.ReferenceFEs: _num_facets
using Gridap.ReferenceFEs: _num_vertices
using Gridap.ReferenceFEs: _num_edges
using Gridap.ReferenceFEs: _get_facedims
using Gridap.ReferenceFEs: _get_offsets
using Gridap.ReferenceFEs: _get_offset
using Gridap.ReferenceFEs: _find_unique_with_indices

using Gridap.Arrays: Reindexed
using Gridap.Arrays: IdentityVector

import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Arrays: reindex
import Gridap.Arrays: get_array
import Gridap.Arrays: lazy_append
import Gridap.CellData: CellField
import Gridap.CellData: CellQuadrature
import Gridap.CellData: QPointCellField
import Gridap.CellData: get_cell_map

import Gridap.Fields: field_cache
import Gridap.Fields: evaluate_field!
import Gridap.Fields: evaluate_field_array
import Gridap.Fields: gradient
import Gridap.Fields: grad2curl
import Gridap.Helpers: operate

import Gridap.Integration: get_coordinates
import Gridap.Integration: get_weights

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
import Gridap.Fields: integrate

import Gridap.Arrays: apply_kernel!

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
export get_normal_vector
export test_triangulation
export restrict
export get_physical_coordinate
export get_cell_id
export cell_measure
export get_cell_map

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
export get_volume_triangulation
export get_face_to_cell
export get_face_to_lface
export get_face_to_cell_map
export get_face_to_face
export get_cell_around
export test_boundary_triangulation

export GenericBoundaryTriangulation

export DiscreteModelPortion

export SkeletonPair
export SkeletonTriangulation
export InterfaceTriangulation
export get_left_boundary
export get_right_boundary

export RestrictedDiscreteModel

export AppendedTriangulation

include("GridTopologies.jl")

include("GridTopologyMocks.jl")

include("UnstructuredGridTopologies.jl")

include("SkeletonPairs.jl")

include("Triangulations.jl")

include("Grids.jl")

include("GridMocks.jl")

include("RestrictedTriangulations.jl")

include("TriangulationPortions.jl")

include("GridPortions.jl")

include("UnstructuredGrids.jl")

include("CartesianGrids.jl")

include("FaceLabelings.jl")

include("DiscreteModels.jl")

include("DiscreteModelPortions.jl")

include("DiscreteModelMocks.jl")

include("UnstructuredDiscreteModels.jl")

include("CartesianDiscreteModels.jl")

include("BoundaryTriangulations.jl")

include("GenericBoundaryTriangulations.jl")

include("SkeletonTriangulations.jl")

include("AppendedTriangulations.jl")

include("RestrictedDiscreteModels.jl")

end # module
