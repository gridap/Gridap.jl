module NodalReferenceFEsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Fields

D = 2
T = Float64
order = 1
prebasis = MonomialBasis{D}(T,order)

polytope = QUAD
x = get_vertex_coordinates(polytope)
dofs = LagrangianDofBasis(T,x)

ndofs = 4

face_own_dofs = [[1],[2],[3],[4],Int[],Int[],Int[],Int[],Int[]]
a = [1]
b = Int[]
face_own_dofs_permutations = [[a],[a],[a],[a],[b,b],[b,b],[b,b],[b,b],fill(b,8)]
face_dofs = [[1],[2],[3],[4],[1,2],[3,4],[1,3],[2,4],[1,2,3,4]]
metadata = nothing

struct MockConformity <: Conformity end

reffe = GenericRefFE(
  ndofs, polytope, prebasis, dofs,
  MockConformity(), metadata, face_dofs)

function ReferenceFEs.get_face_own_dofs(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs
end

function ReferenceFEs.get_face_own_dofs_permutations(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs_permutations
end

node_coordinates = Point{2,Float64}[(0,0),(1,0),(0,1),(1,1)]
node_and_comp_to_dof = [1,2,3,4]
face_own_nodes = face_own_dofs
face_own_nodes_permutations = face_own_dofs_permutations
face_nodes = face_dofs

reffe = GenericNodalCartesianRefFE(
  reffe,
  node_coordinates,
  node_and_comp_to_dof,
  face_own_nodes,
  face_own_nodes_permutations,
  face_nodes)

test_nodal_reference_fe(reffe)

@test get_dof_to_node(reffe) == [1, 2, 3, 4]
@test get_own_nodes_permutations(reffe) == get_face_own_nodes_permutations(reffe)[end]

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_nodes(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_dofs_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]
@test get_face_own_nodes_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]


end # module
