
"""
    abstract type LagrangianRefFE{D} <: ReferenceFE{D}

Abstract type representing a Lagrangian reference FE. Lagrangian in the sense that
`get_dof_basis` returns an instance of `LagrangianDofBasis`.
The interface for this type is defined with the methods of `ReferenceFE`
plus the following ones

- [`get_face_own_nodes(reffe::LagrangianRefFE,conf::Conformity)`](@ref)
- [`get_face_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity)`](@ref)
- [`get_face_nodes(reffe::LagrangianRefFE)`](@ref)

"""
abstract type LagrangianRefFE{D} <: ReferenceFE{D} end

"""
    get_face_own_nodes(reffe::LagrangianRefFE,conf::Conformity)
"""
function get_face_own_nodes(reffe::LagrangianRefFE,conf::Conformity)
  @abstractmethod
end

"""
    get_face_own_nodes(reffe::LagrangianRefFE)
"""
function get_face_own_nodes(reffe::LagrangianRefFE)
  conf = Conformity(reffe)
  get_face_own_nodes(reffe,conf)
end

function get_face_own_nodes(reffe::LagrangianRefFE,conf::L2Conformity)
  _get_face_own_nodes_l2(reffe)
end

function _get_face_own_nodes_l2(reffe::LagrangianRefFE)
  p = get_polytope(reffe)
  r = [Int[] for i in 1:num_faces(p)]
  r[end] = collect(1:num_nodes(reffe))
  r
end

"""
    get_face_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity)
"""
function get_face_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity)
  face_own_nodes = get_face_own_nodes(reffe,conf)
  _trivial_face_own_dofs_permutations(face_own_nodes)
end

"""
    get_face_own_nodes_permutations(reffe::LagrangianRefFE)
"""
function get_face_own_nodes_permutations(reffe::LagrangianRefFE)
  conf = Conformity(reffe)
  get_face_own_nodes_permutations(reffe,conf)
end

"""
    get_face_nodes(reffe::LagrangianRefFE)
"""
function get_face_nodes(reffe::LagrangianRefFE)
  @abstractmethod
end

# Tester

"""
    test_lagrangian_reference_fe(reffe::LagrangianRefFE)
"""
function test_lagrangian_reference_fe(reffe::LagrangianRefFE)
  conf = Conformity(reffe)
  @test isa(conf,Conformity)
  test_lagrangian_reference_fe(reffe,conf)
end

function test_lagrangian_reference_fe(reffe::LagrangianRefFE,conf::Conformity)
  test_reference_fe(reffe,conf)
  D = num_dims(reffe)
  node_coordinates = get_node_coordinates(reffe)
  @test isa(node_coordinates,Vector{<:Point{D}})
  @test length(node_coordinates) == num_nodes(reffe)
  node_and_comp_to_dof = get_node_and_comp_to_dof(reffe)
  @test isa(get_dof_to_node(reffe),Vector{Int})
  @test isa(get_dof_to_comp(reffe),Vector{Int})
  @test isa(node_and_comp_to_dof,Vector)
  dof_to_node = get_dof_to_node(reffe)
  @test isa(dof_to_node,Vector{Int})
  @test isa(get_face_own_nodes(reffe,conf),Vector{Vector{Int}})
  @test isa(get_face_own_nodes_permutations(reffe,conf),Vector{Vector{Vector{Int}}})
  @test isa(get_face_nodes(reffe),Vector{Vector{Int}})
end

# Default API

"""
    get_node_coordinates(reffe::LagrangianRefFE)
"""
function get_node_coordinates(reffe::LagrangianRefFE)
  dofs = get_dof_basis(reffe)
  dofs.nodes
end

"""
    num_nodes(reffe::LagrangianRefFE)
"""
num_nodes(reffe::LagrangianRefFE) = length(get_node_coordinates(reffe))

"""
    get_node_and_comp_to_dof(reffe::LagrangianRefFE)
"""
function get_node_and_comp_to_dof(reffe::LagrangianRefFE)
  dofs = get_dof_basis(reffe)
  dofs.node_and_comp_to_dof
end

"""
    get_dof_to_node(reffe::LagrangianRefFE)
"""
function get_dof_to_node(reffe::LagrangianRefFE)
  dofs = get_dof_basis(reffe)
  dofs.dof_to_node
end

"""
    get_dof_to_comp(reffe::LagrangianRefFE)
"""
function get_dof_to_comp(reffe::LagrangianRefFE)
  dofs = get_dof_basis(reffe)
  dofs.dof_to_comp
end

"""
    get_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity)
"""
function get_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity)
  n = num_faces(reffe)
  get_face_own_nodes_permutations(reffe,conf)[n]
end

"""
    get_own_nodes_permutations(reffe::LagrangianRefFE)
"""
function get_own_nodes_permutations(reffe::LagrangianRefFE)
  conf = Conformity(reffe)
  get_own_nodes_permutations(reffe,conf)
end

"""
    get_vertex_node(reffe::LagrangianRefFE,conf::Conformity) -> Vector{Int}
"""
function get_vertex_node(reffe::LagrangianRefFE,conf::Conformity)
  d = 0
  p = get_polytope(reffe)
  range = get_dimranges(p)[d+1]
  vertex_to_nodes = get_face_own_nodes(reffe,conf)[range]
  map(first, vertex_to_nodes)
end

"""
    get_vertex_node(reffe::LagrangianRefFE) -> Vector{Int}
"""
function get_vertex_node(reffe::LagrangianRefFE)
  conf = Conformity(reffe)
  get_vertex_node(reffe,conf)
end

"""
    get_face_own_nodes(reffe::LagrangianRefFE,conf::Conformity,d::Integer)
"""
function get_face_own_nodes(reffe::LagrangianRefFE,conf::Conformity,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_own_nodes(reffe,conf)[range]
end

"""
    get_face_own_nodes(reffe::LagrangianRefFE,d::Integer)
"""
function get_face_own_nodes(reffe::LagrangianRefFE,d::Integer)
  conf = Conformity(reffe)
  get_face_own_nodes(reffe,conf,d)
end

"""
    get_face_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity,d::Integer)
"""
function get_face_own_nodes_permutations(reffe::LagrangianRefFE,conf::Conformity,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_own_nodes_permutations(reffe,conf)[range]
end

"""
    get_face_own_nodes_permutations(reffe::LagrangianRefFE,d::Integer)
"""
function get_face_own_nodes_permutations(reffe::LagrangianRefFE,d::Integer)
  conf = Conformity(reffe)
  get_face_own_nodes_permutations(reffe,conf,d)
end

"""
    get_face_nodes(reffe::LagrangianRefFE,d::Integer)
"""
function get_face_nodes(reffe::LagrangianRefFE,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_nodes(reffe)[range]
end

function get_face_own_dofs(reffe::LagrangianRefFE,conf::Conformity)
  dofs = get_dof_basis(reffe)
  face_own_nodes = get_face_own_nodes(reffe,conf)
  face_own_dofs = _generate_face_own_dofs(face_own_nodes, dofs.node_and_comp_to_dof)
  face_own_dofs
end

function _generate_face_own_dofs(face_own_nodes, node_and_comp_to_dof)
  faces = 1:length(face_own_nodes)
  T = eltype(node_and_comp_to_dof)
  comps = 1:num_components(T)
  face_own_dofs = [Int[] for i in faces]
  for face in faces
    nodes = face_own_nodes[face]
    # Node major
    for comp in comps
      for node in nodes
        comp_to_dofs = node_and_comp_to_dof[node]
        dof = comp_to_dofs[comp]
        push!(face_own_dofs[face],dof)
      end
    end
  end

  face_own_dofs

end

function get_face_own_dofs_permutations(reffe::LagrangianRefFE,conf::Conformity)
  dofs = get_dof_basis(reffe)
  face_own_nodes_permutations = get_face_own_nodes_permutations(reffe,conf)
  face_own_nodes = get_face_own_nodes(reffe,conf)
  face_own_dofs = get_face_own_dofs(reffe,conf)
  face_own_dofs_permutations = _generate_face_own_dofs_permutations(
    face_own_nodes_permutations, dofs.node_and_comp_to_dof, face_own_nodes, face_own_dofs)
  face_own_dofs_permutations
end

function  _generate_face_own_dofs_permutations(
  face_own_nodes_permutations, node_and_comp_to_dof, face_own_nodes, face_own_dofs)

  T = eltype(node_and_comp_to_dof)
  ncomps = num_components(T)

  face_own_dofs_permutations = Vector{Vector{Int}}[]
  for  (face, pindex_to_inode_to_pinode) in enumerate(face_own_nodes_permutations)
    idof_to_dof = face_own_dofs[face]
    inode_to_node = face_own_nodes[face]
    pindex_to_idof_to_pidof = Vector{Int}[]
    for inode_to_pinode in pindex_to_inode_to_pinode
      ninodes = length(inode_to_pinode)
      nidofs = ncomps*ninodes
      idof_to_pidof = fill(INVALID_PERM,nidofs)
      for (inode,ipnode) in enumerate(inode_to_pinode)
        if ipnode == INVALID_PERM
          continue
        end
        node = inode_to_node[inode]
        pnode = inode_to_node[ipnode]
        comp_to_pdof = node_and_comp_to_dof[pnode]
        comp_to_dof = node_and_comp_to_dof[node]
        for comp in 1:ncomps
          dof = comp_to_dof[comp]
          pdof = comp_to_pdof[comp]
          idof = findfirst(i->i==dof,idof_to_dof)
          ipdof = findfirst(i->i==pdof,idof_to_dof)
          idof_to_pidof[idof] = ipdof
        end
      end
      push!(pindex_to_idof_to_pidof,idof_to_pidof)
    end
    push!(face_own_dofs_permutations,pindex_to_idof_to_pidof)
  end
  face_own_dofs_permutations
end

function get_face_own_dofs(reffe::LagrangianRefFE,conf::L2Conformity)
  _get_face_own_dofs_l2(reffe)
end

function get_face_own_dofs_permutations(reffe::LagrangianRefFE,conf::L2Conformity)
  face_own_dofs = get_face_own_dofs(reffe,conf)
  _trivial_face_own_dofs_permutations(face_own_dofs)
end

"""
    (==)(a::LagrangianRefFE{D}, b::LagrangianRefFE{D}) where D
"""
function (==)(a::LagrangianRefFE{D}, b::LagrangianRefFE{D}) where D
  t = true
  pa = get_polytope(a)
  pb = get_polytope(b)
  t = t && (pa == pb)
  xa = get_node_coordinates(a)
  xb = get_node_coordinates(b)
  t = t && (xa ≈ xb)
  expsa = get_exponents(get_prebasis(a))
  expsb = get_exponents(get_prebasis(b))
  t = t && (expsa == expsb)
  facedofsa = get_face_dofs(a)
  facedofsb = get_face_dofs(b)
  t = t && (facedofsa == facedofsb)
  ia = get_node_and_comp_to_dof(a)
  ib = get_node_and_comp_to_dof(b)
  t = t && (ia == ib)
  t
end

function (==)(a::LagrangianRefFE, b::LagrangianRefFE)
  false
end

"""
    get_order(reffe::LagrangianRefFE)
"""
function get_order(reffe::LagrangianRefFE)
  get_order(get_prebasis(reffe))
end

"""
    get_orders(reffe::LagrangianRefFE)
"""
function get_orders(reffe::LagrangianRefFE)
  get_orders(get_prebasis(reffe))
end

# Generic implementation
"""
  struct GenericLagrangianRefFE{C,D} <: LagrangianRefFE{D}
    reffe::GenericRefFE{C,D}
    face_nodes::Vector{Vector{Int}}
  end
"""
struct GenericLagrangianRefFE{C,D} <: LagrangianRefFE{D}
  reffe::GenericRefFE{C,D}
  face_nodes::Vector{Vector{Int}}
end

# LagrangianRefFE

get_face_nodes(reffe::GenericLagrangianRefFE) = reffe.face_nodes

# Reffe

num_dofs(reffe::GenericLagrangianRefFE) = num_dofs(reffe.reffe)

get_polytope(reffe::GenericLagrangianRefFE) = get_polytope(reffe.reffe)

get_prebasis(reffe::GenericLagrangianRefFE) = get_prebasis(reffe.reffe)

get_dof_basis(reffe::GenericLagrangianRefFE) = get_dof_basis(reffe.reffe)

Conformity(reffe::GenericLagrangianRefFE) = Conformity(reffe.reffe)

get_face_dofs(reffe::GenericLagrangianRefFE) = get_face_dofs(reffe.reffe)

get_shapefuns(reffe::GenericLagrangianRefFE) = get_shapefuns(reffe.reffe)

