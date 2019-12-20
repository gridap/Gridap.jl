
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


# Dafault API

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
end


# Generic implementation

"""
  struct GenericNodalRefFE{D,T,V} <: NodalReferenceFE{D}
    reffe::GenericRefFE{D}
    node_coordinates::Vector{Point{D,T}}
    node_and_comp_to_dof::Vector{V}
  end
"""
struct GenericNodalRefFE{D,T,V} <: NodalReferenceFE{D}
  reffe::GenericRefFE{D}
  node_coordinates::Vector{Point{D,T}}
  node_and_comp_to_dof::Vector{V}
end

get_node_coordinates(reffe::GenericNodalRefFE) = reffe.node_coordinates

get_node_and_comp_to_dof(reffe::GenericNodalRefFE) = reffe.node_and_comp_to_dof

num_dofs(reffe::GenericNodalRefFE) = reffe.reffe.ndofs

get_polytope(reffe::GenericNodalRefFE) = reffe.reffe.polytope

get_prebasis(reffe::GenericNodalRefFE) = reffe.reffe.prebasis

get_dof_basis(reffe::GenericNodalRefFE) = reffe.reffe.dofs

get_face_own_dofs(reffe::GenericNodalRefFE) = reffe.reffe.face_own_dofs

get_face_own_dofs_permutations(reffe::GenericNodalRefFE) = reffe.reffe.face_own_dofs_permutations

get_face_dofs(reffe::GenericNodalRefFE) = reffe.reffe.face_own_dofs

get_shapefuns(reffe::GenericNodalRefFE) = reffe.reffe.shapefuns

