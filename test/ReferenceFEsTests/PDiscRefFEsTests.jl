module PDiscRefFEsTests

using Test
using Gridap.TensorValues
using Gridap.ReferenceFEs

T = Float64
reffe = LagrangianRefFE(T,QUAD,2,space=:P)
@test reffe == ReferenceFE(QUAD,:S,2,2) # r=2,k=D=2

T = Float64
reffe = LagrangianRefFE(T,HEX,2,space=:P)
@test reffe == ReferenceFE(HEX,:S,2,3,T) # r=2,k=D=3

T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,QUAD,2,space=:P)
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
