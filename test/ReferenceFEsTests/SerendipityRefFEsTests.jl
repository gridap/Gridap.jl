module SerendipityRefFEsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Io

nodal = true

reffe = SerendipityRefFE(Float64,QUAD,2)
@test reffe == ReferenceFE(QUAD,:S,2,0; nodal)
@test get_node_coordinates(reffe) == Point{2,Float64}[
  (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0),
  (0.5, 0.0), (0.5, 1.0), (0.0, 0.5), (1.0, 0.5)]
@test get_face_own_nodes(reffe) == [[1], [2], [3], [4], [5], [6], [7], [8], Int[]]
@test get_face_own_dofs(reffe) == Vector{Int}[[1], [2], [3], [4], [5], [6], [7], [8],[]]
@test get_face_dofs(reffe) == Vector{Int}[
  [1], [2], [3], [4],
  [1,2,5], [3,4,6], [1,3,7], [2,4,8],
  [1,2,3,4,5,6,7,8]
]
test_lagrangian_reference_fe(reffe)

reffer = ReferenceFE(QUAD,serendipity,2)
@test typeof(reffe) == typeof(reffer)
@test get_node_coordinates(reffe) == get_node_coordinates(reffer)
@test get_face_own_nodes(reffe) == get_face_own_nodes(reffer)

reffe = SerendipityRefFE(Float64,QUAD,4)
@test reffe == ReferenceFE(QUAD,:S,4,0; nodal)
@test get_face_own_nodes(reffe) == [[1], [2], [3], [4], [5, 6, 7], [8, 9, 10], [11, 12, 13], [14, 15, 16], [17]]
test_lagrangian_reference_fe(reffe)

reffe = SerendipityRefFE(Float64,HEX,2)
@test reffe == ReferenceFE(HEX,:S,2,0; nodal)
@test get_face_own_nodes(reffe) == [
  [1], [2], [3], [4], [5], [6], [7], [8],
  [9], [10], [11], [12], [13], [14], [15], [16],
  [17], [18], [19], [20], Int[], Int[], Int[], Int[], Int[], Int[], Int[]]
test_lagrangian_reference_fe(reffe)

p = get_polytope(reffe)
test_polytope(p)
@test LagrangianRefFE(p) == HEX8

poly_type=Polynomials.ModalC0

order = 4
reffe = SerendipityRefFE(Float64,HEX,order)
test_lagrangian_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; nodal)
@test get_order(reffe) == order
@test get_orders(reffe) == (order,order,order)
@test is_P(reffe) == false
@test is_Q(reffe) == false
@test is_S(reffe) == true

@test_throws "hierarchical" reffe = SerendipityRefFE(Float64,HEX,order; poly_type=Bernstein)

reffe = SerendipityRefFE(Float64,HEX,order; poly_type)
test_lagrangian_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, nodal)

space=:S

reffe = LagrangianRefFE(Float64,HEX,order; space)
test_lagrangian_reference_fe(reffe)
@test get_order(reffe) == order
@test get_orders(reffe) == (order,order,order)
@test is_P(reffe) == false
@test is_Q(reffe) == false
@test is_S(reffe) == true

reffe = LagrangianRefFE(Float64,HEX,order; space, poly_type)
test_lagrangian_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, nodal)

reffe = SerendipityRefFE(Float64,QUAD,(3,3))
test_lagrangian_reference_fe(reffe)
@test reffe == ReferenceFE(QUAD,:S,3,0,Float64; nodal)
@test reffe == from_dict(LagrangianRefFE,to_dict(reffe))

end # module
