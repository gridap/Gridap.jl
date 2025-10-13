module PDiscRefFEsTests

using Test
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs

T = Float64
reffe = LagrangianRefFE(T,QUAD,3,space=:P)
@test reffe == ReferenceFE(QUAD,:S,3,2) # r=2,k=D=2

reffe = LagrangianRefFE(T,QUAD,3,space=:P; poly_type=Bernstein)
@test reffe == ReferenceFE(QUAD,:S,3,2; poly_type=Bernstein) # r=2,k=D=2

reffe = LagrangianRefFE(T,QUAD,3,space=:P; poly_type=Bernstein, sh_is_pb=true)
@test reffe == ReferenceFE(QUAD,:S,3,2; poly_type=Bernstein, sh_is_pb=true) # r=2,k=D=2

reffe = LagrangianRefFE(T,QUAD,3,space=:P; sh_is_pb=true)
@test reffe == ReferenceFE(QUAD,:S,3,2; sh_is_pb=true) # r=2,k=D=2

T = Float64
reffe = LagrangianRefFE(T,HEX,2,space=:P)
@test reffe == ReferenceFE(HEX,:S,2,3,T) # r=2,k=D=3

reffe = LagrangianRefFE(T,HEX,2,space=:P; poly_type=Bernstein)
@test reffe == ReferenceFE(HEX,:S,2,3,T; poly_type=Bernstein) # r=2,k=D=3

reffe = LagrangianRefFE(T,HEX,2,space=:P; poly_type=Bernstein, sh_is_pb=true)
@test reffe == ReferenceFE(HEX,:S,2,3,T; poly_type=Bernstein, sh_is_pb=true) # r=2,k=D=3

reffe = LagrangianRefFE(T,HEX,2,space=:P; sh_is_pb=true)
@test reffe == ReferenceFE(HEX,:S,2,3,T; sh_is_pb=true) # r=2,k=D=3

T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,QUAD,2,space=:P)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

reffe = LagrangianRefFE(T,QUAD,2,space=:P, poly_type=Bernstein,)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

reffe = LagrangianRefFE(T,QUAD,2,space=:P, poly_type=Bernstein, sh_is_pb=true)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

reffe = LagrangianRefFE(T,QUAD,2,space=:P, sh_is_pb=true)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,HEX,3,space=:P)
@test Conformity(reffe) == L2Conformity()
test_lagrangian_reference_fe(reffe)
@test is_n_cube(get_polytope(reffe))

#using Gridap.Visualization
#writevtk(reffe,"reffe")

end # module
