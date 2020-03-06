
"""
    abstract type NodalReferenceFE{D} <: ReferenceFE{D}

Abstract type representing a node-based reference FE.
We understand a node-based reference FE as one that
uses the concept of node to locate dofs on the underlying polytope.
Here, nodal-based does not necessary mean an interpolatory reference FE.
We only assume that each dof is assigned to a node, whereas
several dofs can share a same node in general.

The interface for this type is defined with the methods of `ReferenceFE`
plus the following ones

- [`get_node_coordinates(reffe::NodalReferenceFE)`](@ref)
- [`get_node_and_comp_to_dof(reffe::NodalReferenceFE)`](@ref)
- [`get_face_own_nodes(reffe::NodalReferenceFE)`](@ref)
- [`get_face_own_nodes_permutations(reffe::NodalReferenceFE)`](@ref)
- [`get_face_nodes(reffe::NodalReferenceFE)`](@ref)

"""
abstract type NodalReferenceFE{D} <: ReferenceFE{D} end

"""
    get_node_coordinates(reffe::NodalReferenceFE)
"""
function get_node_coordinates(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_node_and_comp_to_dof(reffe::NodalReferenceFE)
"""
function get_node_and_comp_to_dof(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_face_own_nodes(reffe::NodalReferenceFE)
"""
function get_face_own_nodes(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_face_own_nodes_permutations(reffe::NodalReferenceFE)
"""
function get_face_own_nodes_permutations(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_face_nodes(reffe::NodalReferenceFE)
"""
function get_face_nodes(reffe::NodalReferenceFE)
  @abstractmethod
end

# Tester

"""
    test_nodal_reference_fe(reffe::NodalReferenceFE)
"""
function test_nodal_reference_fe(reffe::NodalReferenceFE)
  test_reference_fe(reffe)
  D = num_dims(reffe)
  node_coordinates = get_node_coordinates(reffe)
  @test isa(node_coordinates,Vector{<:Point{D}})
  @test length(node_coordinates) == num_nodes(reffe)
  node_and_comp_to_dof = get_node_and_comp_to_dof(reffe)
  @test isa(node_and_comp_to_dof,Vector)
  dof_to_node = get_dof_to_node(reffe)
  @test isa(dof_to_node,Vector{Int})
  @test isa(get_face_own_nodes(reffe),Vector{Vector{Int}})
  @test isa(get_face_own_nodes_permutations(reffe),Vector{Vector{Vector{Int}}})
  @test isa(get_face_nodes(reffe),Vector{Vector{Int}})
end

# Default API

"""
    num_nodes(reffe::NodalReferenceFE)
"""
num_nodes(reffe::NodalReferenceFE) = length(get_node_coordinates(reffe))

"""
    get_dof_to_node(reffe::NodalReferenceFE)
"""
function get_dof_to_node(reffe::NodalReferenceFE)
  ndofs = num_dofs(reffe)
  dof_to_node = zeros(Int,ndofs)
  nnodes = num_nodes(reffe)
  node_and_comp_to_dof = get_node_and_comp_to_dof(reffe)
  for node in 1:nnodes
    comp_to_dof = node_and_comp_to_dof[node]
    for dof in comp_to_dof
      dof_to_node[dof] = node
    end
  end
  dof_to_node
end

"""
    get_own_nodes_permutations(reffe::NodalReferenceFE)
"""
function get_own_nodes_permutations(reffe::NodalReferenceFE)
  n = num_faces(reffe)
  get_face_own_nodes_permutations(reffe)[n]
end

"""
    get_vertex_node(reffe::NodalReferenceFE) -> Vector{Int}
"""
function get_vertex_node(reffe::NodalReferenceFE)
  d = 0
  p = get_polytope(reffe)
  range = get_dimranges(p)[d+1]
  vertex_to_nodes = get_face_own_nodes(reffe)[range]
  map(first, vertex_to_nodes)
end

"""
    get_face_own_nodes(reffe::NodalReferenceFE,d::Integer)
"""
function get_face_own_nodes(reffe::NodalReferenceFE,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_own_nodes(reffe)[range]
end

"""
    get_face_own_nodes_permutations(reffe::NodalReferenceFE,d::Integer)
"""
function get_face_own_nodes_permutations(reffe::NodalReferenceFE,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_own_nodes_permutations(reffe)[range]
end

"""
    get_face_nodes(reffe::NodalReferenceFE,d::Integer)
"""
function get_face_nodes(reffe::NodalReferenceFE,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_nodes(reffe)[range]
end

# Generic implementation

"""
  struct GenericNodalCartesianRefFE{D,T,V} <: NodalReferenceFE{D}
    reffe::GenericRefFE{D}
    node_coordinates::Vector{Point{D,T}}
    node_and_comp_to_dof::Vector{V}
    face_own_nodes::Vector{Vector{Int}}
    face_own_nodes_permutations::Vector{Vector{Vector{Int}}}
    face_nodes::Vector{Vector{Int}}
  end
"""
struct GenericNodalCartesianRefFE{D,T,V} <: NodalReferenceFE{D}
  reffe::GenericRefFE{D}
  node_coordinates::Vector{Point{D,T}}
  node_and_comp_to_dof::Vector{V}
  face_own_nodes::Vector{Vector{Int}}
  face_own_nodes_permutations::Vector{Vector{Vector{Int}}}
  face_nodes::Vector{Vector{Int}}
end

# NodalReffe

get_node_coordinates(reffe::GenericNodalCartesianRefFE) = reffe.node_coordinates

get_node_and_comp_to_dof(reffe::GenericNodalCartesianRefFE) = reffe.node_and_comp_to_dof

get_face_own_nodes(reffe::GenericNodalCartesianRefFE) = reffe.face_own_nodes

get_face_own_nodes_permutations(reffe::GenericNodalCartesianRefFE) = reffe.face_own_nodes_permutations

get_face_nodes(reffe::GenericNodalCartesianRefFE) = reffe.face_nodes

# Reffe

num_dofs(reffe::GenericNodalCartesianRefFE) = reffe.reffe.ndofs

get_polytope(reffe::GenericNodalCartesianRefFE) = reffe.reffe.polytope

get_prebasis(reffe::GenericNodalCartesianRefFE) = reffe.reffe.prebasis

get_dof_basis(reffe::GenericNodalCartesianRefFE) = reffe.reffe.dofs

get_face_own_dofs(reffe::GenericNodalCartesianRefFE) = reffe.reffe.face_own_dofs

get_face_own_dofs_permutations(reffe::GenericNodalCartesianRefFE) = reffe.reffe.face_own_dofs_permutations

get_face_dofs(reffe::GenericNodalCartesianRefFE) = reffe.reffe.face_own_dofs

get_shapefuns(reffe::GenericNodalCartesianRefFE) = reffe.reffe.shapefuns
