module GeometricLagrangianRefFEsTests

using Test
using Gridap.Io
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using JSON

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
extrusion = Tuple(QUAD.extrusion)

dofs = LagrangianDofBasis(VectorValue{3,Float64},TET,1)
@test dofs.nodes == Point{3,Float64}[(0,0,0), (1,0,0), (0,1,0), (0,0,1)]
@test dofs.node_and_comp_to_dof == VectorValue{3,Int}[(1,5,9), (2,6,10), (3,7,11), (4,8,12)]

dofs = LagrangianDofBasis(Float64,WEDGE,(2,2,2))
r = Point{3,Float64}[
  (0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),
  (0.0,0.0,1.0),(1.0,0.0,1.0),(0.0,1.0,1.0),
  (0.5,0.0,0.0),(0.5,0.0,1.0),(0.0,0.5,0.0),
  (0.5,0.5,0.0),(0.0,0.5,1.0),(0.5,0.5,1.0),
  (0.0,0.0,0.5),(1.0,0.0,0.5),(0.0,1.0,0.5),
  (0.5,0.0,0.5),(0.0,0.5,0.5),(0.5,0.5,0.5)]
@test dofs.nodes == r

dofs = LagrangianDofBasis(VectorValue{2,Int},VERTEX,())
@test dofs.node_and_comp_to_dof == VectorValue{2,Int}[(1,2)]

b = MonomialBasis(VectorValue{2,Int},VERTEX,())
@test length(b) == 2
@test evaluate(b,Point{0,Int}[(),()]) == VectorValue{2,Int}[(1, 0) (0, 1); (1, 0) (0, 1)]

reffe = LagrangianRefFE(VectorValue{2,Int},VERTEX,())
@test get_face_own_nodes(reffe) == [[1]]
@test get_face_own_dofs(reffe) == [[1,2]]
@test get_face_own_dofs_permutations(reffe) == [[[1, 2]]]
test_lagrangian_reference_fe(reffe)
@test ReferenceFE{0}(reffe,1) === reffe

reffe = LagrangianRefFE(VectorValue{2,Float64},SEGMENT,(2,))
@test get_face_own_dofs(reffe) == [[1, 4], [2, 5], [3, 6]]
test_lagrangian_reference_fe(reffe)

reffe = LagrangianRefFE(VectorValue{2,Float64},TRI,3)
@test get_face_own_dofs(reffe) == [[1, 11], [2, 12], [3, 13], [4, 5, 14, 15], [6, 7, 16, 17], [8, 9, 18, 19], [10, 20]]
test_lagrangian_reference_fe(reffe)

reffe = ReferenceFE(TRI,lagrangian,VectorValue{2,Float64},3)
@test get_face_own_dofs(reffe) == [[1, 11], [2, 12], [3, 13], [4, 5, 14, 15], [6, 7, 16, 17], [8, 9, 18, 19], [10, 20]]
test_lagrangian_reference_fe(reffe)

@test Lagrangian() == lagrangian

reffe = LagrangianRefFE(Float64,HEX,2)
test_lagrangian_reference_fe(reffe)

reffe = LagrangianRefFE(Float64,WEDGE,(1,1,2))
test_lagrangian_reference_fe(reffe)
refface = ReferenceFE{1}(reffe,3)
@test get_face_own_dofs(refface) == [[1], [2], []]
refface = ReferenceFE{1}(reffe,8)
@test get_face_own_dofs(refface) == [[1], [2], [3]]

orders = (4,)
reffe = LagrangianRefFE(VectorValue{2,Float64},SEGMENT,orders)
@test get_own_nodes_permutations(reffe) == [[1, 2, 3], [3, 2, 1]]
@test get_own_dofs_permutations(reffe) == [[1, 2, 3, 4, 5, 6], [3, 2, 1, 6, 5, 4]]

orders = (2,3)
reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,orders)
@test get_own_nodes_permutations(reffe) ==[[1, 2], [0, 0], [1, 2], [0, 0], [0, 0], [2, 1], [0, 0], [2, 1]]
@test get_own_dofs_permutations(reffe) == [
  [1, 2, 3, 4], [0, 0, 0, 0], [1, 2, 3, 4], [0, 0, 0, 0],
  [0, 0, 0, 0], [2, 1, 4, 3], [0, 0, 0, 0], [2, 1, 4, 3]]

reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,2)
@test get_node_and_comp_to_dof(reffe) == VectorValue{2,Int}[
  (1, 10), (2, 11), (3, 12), (4, 13), (5, 14), (6, 15), (7, 16), (8, 17), (9, 18)]

reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,2)
d = 1
@test get_face_own_dofs(reffe,d) == [[5, 14], [6, 15], [7, 16], [8, 17]]
@test get_face_own_nodes(reffe,d) == [[5], [6], [7], [8]]

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
@test num_nodes(reffe) == 1
test_lagrangian_reference_fe(reffe)
@test Conformity(reffe) == L2Conformity()

@test get_face_own_nodes(reffe) == Vector{Int}[[], [], [], [], [], [], [], [], [1]]
#@test get_face_own_nodes_permutations(reffe) == Vector{Vector{Int}}[
#  [[]],[[]],[[]],[[]],[[],[]],[[],[]],[[],[]],[[],[]],[[1],[1],[1],[1],[1],[1],[1],[1]] ]
#@test get_own_nodes_permutations(reffe) == [[1], [1], [1], [1], [1], [1], [1], [1]]
@test get_face_own_nodes_permutations(reffe) == Vector{Vector{Int}}[
  [[]],[[]],[[]],[[]],[[]],[[]],[[]],[[]],[[1]] ]
@test get_own_nodes_permutations(reffe) == [[1]]

reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,orders)
@test num_dofs(reffe) == 2
@test num_nodes(reffe) == 1
test_lagrangian_reference_fe(reffe)

@test get_face_own_dofs(reffe) == Array{Int,1}[[], [], [], [], [], [], [], [], [1, 2]]
#@test get_face_own_dofs_permutations(reffe) == [
#  [[]],[[]],[[]],[[]],[[],[]],[[],[]],[[],[]],[[],[]],[[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2]]]
@test get_face_own_dofs_permutations(reffe) == [
  [[]],[[]],[[]],[[]],[[]],[[]],[[]],[[]],[[1,2]]]

own_nodes = compute_own_nodes(TRI,orders)
@test own_nodes == Point{2,Float64}[(1.0/3,1.0/3)]

# More API

reffe = LagrangianRefFE(Float64,WEDGE,1)

reffes = get_reffaces(ReferenceFE{2}, reffe)
iface_to_ftype = get_face_type(reffe,2)

@test length(reffes) == 2
@test iface_to_ftype == [1, 1, 2, 2, 2]

order = 4
reffe = LagrangianRefFE(Float64,TET,order)
@test get_order(reffe) == order
@test get_orders(reffe) == (order,order,order)
@test is_P(reffe) == true
@test is_Q(reffe) == false
@test is_S(reffe) == false

orders = (1,2)
reffe = LagrangianRefFE(Float64,QUAD,orders)
@test get_order(reffe) == 2
@test get_orders(reffe) == orders
@test is_P(reffe) == false
@test is_Q(reffe) == true
@test is_S(reffe) == false

reffe = LagrangianRefFE(Float64,QUAD,(2,3))
@test reffe == from_dict(LagrangianRefFE,to_dict(reffe))

reffe = LagrangianRefFE(Float64,QUAD,(2,3))

s = to_json(reffe)
@test reffe == from_json(LagrangianRefFE,s)

s = JSON.json(reffe)
@test reffe == from_json(LagrangianRefFE,s)

d = mktempdir()
f = joinpath(d,"reffe.jld2")

to_jld2_file(reffe,f)
@test reffe == from_jld2_file(typeof(reffe),f)

rm(d,recursive=true)

end # module
