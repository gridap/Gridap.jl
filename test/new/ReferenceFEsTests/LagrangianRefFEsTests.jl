module LagrangianRefFEsTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs

orders = (2,3)
b = MonomialBasis(Float64,QUAD,orders)
r = [(0,0), (1,0), (2,0), (0,1), (1,1), (2,1), (0,2), (1,2), (2,2), (0,3), (1,3), (2,3)]
@test get_exponents(b) == r

orders = (1,1,2)
b = MonomialBasis(Float64,WEDGE,orders)
r = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,0,1), (0,1,1), (0,0,2), (1,0,2), (0,1,2)]
@test get_exponents(b) == r

orders = (1,1,1)
b = MonomialBasis(Float64,PYRAMID,orders)
r = [(0,0,0), (1,0,0), (0,1,0), (1,1, 0), (0,0,1)]
@test get_exponents(b) == r

orders = (1,1,1)
b = MonomialBasis(Float64,TET,orders)
r = [(0,0,0), (1,0,0), (0,1,0), (0,0,1)]
@test get_exponents(b) == r

orders = (2,2)
extrusion = Tuple(QUAD.extrusion.array)

dofs = LagrangianDofBasis(VectorValue{3,Float64},TET,1)
@test dofs.nodes == Point{3,Float64}[(0,0,0), (1,0,0), (0,1,0), (0,0,1)]
@test dofs.node_and_comp_to_dof == VectorValue{3,Int}[(1,5,9), (2,6,10), (3,7,11), (4,8,12)]

dofs = LagrangianDofBasis(Float64,WEDGE,(2,2,2))
@test dofs.nodes == Point{3,Float64}[
  (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
  (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0),
  (0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (0.0, 1.0, 0.5),
  (0.5, 0.0, 0.0), (0.5, 0.0, 1.0), (0.0, 0.5, 0.0),
  (0.5, 0.5, 0.0), (0.0, 0.5, 1.0), (0.5, 0.5, 1.0),
  (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.5, 0.5, 0.5)]

dofs = LagrangianDofBasis(VectorValue{2,Int},VERTEX,())
@test dofs.node_and_comp_to_dof == VectorValue{2,Int}[(1,2)]

b = MonomialBasis(VectorValue{2,Int},VERTEX,())
@test evaluate(b,Point{0,Int}[(),()]) == VectorValue{2,Int}[(1, 0) (0, 1); (1, 0) (0, 1)]

reffe = LagrangianRefFE(VectorValue{2,Int},VERTEX,())
@test reffe.face_own_nodeids == [[1]]
@test reffe.data.face_own_dofids == [[1,2]]
test_nodal_reference_fe(reffe,optional=true)
@test ReferenceFE{0}(reffe,1) === reffe

reffe = LagrangianRefFE(VectorValue{2,Float64},SEGMENT,(2,))
@test get_face_own_dofids(reffe) == [[1, 4], [2, 5], [3, 6]]
test_nodal_reference_fe(reffe,optional=true)

reffe = LagrangianRefFE(VectorValue{2,Float64},TRI,3)
@test get_face_own_dofids(reffe) == [[1, 11], [2, 12], [3, 13], [4, 5, 14, 15], [6, 7, 16, 17], [8, 9, 18, 19], [10, 20]]
test_nodal_reference_fe(reffe,optional=true)

reffe = LagrangianRefFE(Float64,HEX,2)
test_nodal_reference_fe(reffe,optional=true)

reffe = LagrangianRefFE(Float64,WEDGE,(1,1,2))
test_nodal_reference_fe(reffe,optional=true)
refface = ReferenceFE{1}(reffe,3)
@test get_face_own_dofids(refface) == [[1], [2], [3]]
refface = ReferenceFE{1}(reffe,4)
@test get_face_own_dofids(refface) == [[1], [2], []]

orders = (4,)
reffe = LagrangianRefFE(VectorValue{2,Float64},SEGMENT,orders)
@test reffe.own_nodes_permutations == [[1, 2, 3], [3, 2, 1]]
@test get_own_dofs_permutations(reffe) == [[1, 2, 3, 4, 5, 6], [3, 2, 1, 6, 5, 4]] 

orders = (2,3)
reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,orders)
@test reffe.own_nodes_permutations ==[[1, 2], [0, 0], [1, 2], [0, 0], [0, 0], [2, 1], [0, 0], [2, 1]] 
@test get_own_dofs_permutations(reffe) == [
  [1, 2, 3, 4], [0, 0, 0, 0], [1, 2, 3, 4], [0, 0, 0, 0],
  [0, 0, 0, 0], [2, 1, 4, 3], [0, 0, 0, 0], [2, 1, 4, 3]]


# 0-order degenerated case

orders = (0,0)

b = compute_monomial_basis(Float64,QUAD,orders)
@test length(b.terms) == 1

own_nodes = compute_own_nodes(QUAD,orders)
@test own_nodes == Point{2,Float64}[(0.5, 0.5)]

face_orders = compute_face_orders(QUAD,SEGMENT,1,orders)
@test face_orders == (0,)

nodes, facenodeids = compute_nodes(QUAD,orders)
@test nodes == Point{2,Float64}[(0.5, 0.5)]
@test facenodeids == [Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[], [1]]

nodeperms = compute_own_nodes_permutations(QUAD,own_nodes)
@test nodeperms == [[1], [1], [1], [1], [1], [1], [1], [1]]

reffaces = compute_lagrangian_reffaces(Float64,QUAD,orders)

reffe = LagrangianRefFE(Float64,QUAD,orders)
@test num_dofs(reffe) == 1

reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,orders)
@test num_dofs(reffe) == 2

own_nodes = compute_own_nodes(TRI,orders)
@test own_nodes == Point{2,Float64}[(1.0/3,1.0/3)]

reffe = LagrangianRefFE(Float64,WEDGE,1)

reffes = get_reffes(ReferenceFE{2}, reffe)

iface_to_ftype = get_face_types(ReferenceFE{2}, reffe)

@test length(reffes) == 2
@test iface_to_ftype == [1, 1, 1, 2, 2]

end # module
