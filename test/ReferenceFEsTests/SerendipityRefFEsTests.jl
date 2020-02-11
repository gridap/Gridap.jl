module SerendipityRefFEsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Io

reffe = SerendipityRefFE(Float64,QUAD,2)
@test reffe.data.dofs.nodes == Point{2,Float64}[
  (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0),
  (0.5, 0.0), (0.5, 1.0), (0.0, 0.5), (1.0, 0.5)]
@test reffe.face_own_nodes == [[1], [2], [3], [4], [5], [6], [7], [8], Int[]]
test_nodal_reference_fe(reffe)

reffe = SerendipityRefFE(Float64,QUAD,4)
@test reffe.face_own_nodes == [[1], [2], [3], [4], [5, 6, 7], [8, 9, 10], [11, 12, 13], [14, 15, 16], [17]] 
test_nodal_reference_fe(reffe)

reffe = SerendipityRefFE(Float64,HEX,2)
@test reffe.face_own_nodes == [
  [1], [2], [3], [4], [5], [6], [7], [8],
  [9], [10], [11], [12], [13], [14], [15], [16],
  [17], [18], [19], [20], Int[], Int[], Int[], Int[], Int[], Int[], Int[]]
test_nodal_reference_fe(reffe)

p = get_polytope(reffe)
test_polytope(p)
@test NodalReferenceFE(p) == HEX8

order = 4
reffe = SerendipityRefFE(Float64,HEX,order)
@test get_order(reffe) == order
@test get_orders(reffe) == (order,order,order)
@test is_P(reffe) == false
@test is_Q(reffe) == false
@test is_S(reffe) == true

reffe = SerendipityRefFE(Float64,QUAD,(3,3))
@test reffe == from_dict(LagrangianRefFE,to_dict(reffe))

end # module
