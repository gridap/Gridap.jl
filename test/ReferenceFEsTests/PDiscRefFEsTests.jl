module PDiscRefFEsTests

using Test
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs

nodal = true

T = Float64
reffe = LagrangianRefFE(T,QUAD,1,space=:P)
@test reffe == ReferenceFE(QUAD,:S,1,2; nodal) # r=1,k=D=2

face_own_dofs = Vector{Int}[[],[],[],[],[],[],[],[],[1,2,3]]
face_dofs = Vector{Int}[[],[],[],[],[],[],[],[],[1,2,3]]
@test get_face_own_dofs(reffe) == face_own_dofs
@test get_face_dofs(reffe) == face_dofs

reffe = LagrangianRefFE(T,QUAD,3,space=:P)
@test reffe == ReferenceFE(QUAD,:S,3,2; nodal) # r=3,k=D=2

reffe = LagrangianRefFE(T,QUAD,3,space=:P; poly_type=Bernstein)
@test reffe == ReferenceFE(QUAD,:S,3,2; poly_type=Bernstein, nodal) # r=3,k=D=2


T = Float64
reffe = LagrangianRefFE(T,HEX,2,space=:P)
@test reffe == ReferenceFE(HEX,:S,2,3,T; nodal) # r=2,k=D=3

reffe = LagrangianRefFE(T,HEX,2,space=:P; poly_type=Bernstein)
@test reffe == ReferenceFE(HEX,:S,2,3,T; poly_type=Bernstein, nodal) # r=2,k=D=3


T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,QUAD,2,space=:P)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

reffe = LagrangianRefFE(T,QUAD,2,space=:P, poly_type=Bernstein,)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))


T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,HEX,3,space=:P)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

end # module
