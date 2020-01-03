"""

The exported names are
$(EXPORTS)
"""
module ReferenceFEs

using Test
using DocStringExtensions
using LinearAlgebra
using Combinatorics

using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Polynomials

import Gridap.Arrays: kernel_cache
import Gridap.Arrays: apply_kernel!
import Gridap.Arrays: kernel_return_type
import Gridap.Fields: evaluate
import Gridap.Polynomials: MonomialBasis

import Gridap.Polynomials: get_order
import Gridap.Polynomials: get_orders

import Gridap.Io: to_dict
import Gridap.Io: from_dict

import Base: ==

export Polytope
export ExtrusionPolytope
export get_extrusion
export get_faces
export get_dimranges
export get_dimrange
export get_vertex_coordinates
export get_facet_normals
export get_facet_orientations
export get_edge_tangents
export get_vertex_permutations
export is_simplex
export is_n_cube
export simplexify
export num_dims
export num_cell_dims
export num_point_dims
export num_vertices
export num_faces
export num_facets
export num_edges
export get_facedims
export get_offsets
export get_offset
export get_face_vertices
export get_reffaces
export get_face_type
export get_bounding_box
export get_face_vertex_permutations
export test_polytope
export VERTEX
export SEGMENT
export TRI
export QUAD
export TET
export HEX
export WEDGE
export PYRAMID
export HEX_AXIS
export TET_AXIS
export INVALID_PERM

export Dof
export evaluate_dof!
export evaluate_dof
export dof_cache
export dof_return_type
export test_dof
export evaluate_dof_array

export ReferenceFE
export GenericRefFE
export get_polytope
export get_prebasis
export get_dof_basis
export get_face_own_dofs
export get_face_own_dofs_permutations
export get_face_dofs
export get_own_dofs_permutations
export get_shapefuns
export compute_shapefuns
export test_reference_fe
export num_dofs

export NodalReferenceFE
export GenericNodalRefFE
export get_face_own_nodes
export get_face_nodes
export get_face_own_nodes_permutations
export get_own_nodes_permutations
export get_node_coordinates
export get_dof_to_node
export get_dof_to_comp
export get_node_and_comp_to_dof
export get_vertex_node
export num_nodes
export test_nodal_reference_fe

export LagrangianRefFE
export LagrangianDofBasis
export compute_monomial_basis
export compute_own_nodes
export compute_face_orders
export compute_nodes
export compute_own_nodes_permutations
export compute_lagrangian_reffaces
export is_first_order
export is_affine
export is_Q
export is_P
export is_S

export VERTEX1
export SEG2
export TRI3
export QUAD4
export TET4
export HEX8

export SerendipityRefFE

include("Polytopes.jl")

include("ExtrusionPolytopes.jl")

include("Dofs.jl")

include("MockDofs.jl")

include("LagrangianDofBases.jl")

include("ReferenceFEInterfaces.jl")

include("NodalReferenceFEs.jl")

include("LagrangianRefFEs.jl")

include("SerendipityRefFEs.jl")

end # module
