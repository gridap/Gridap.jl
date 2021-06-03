module CDLagrangianRefFEsTests

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Test

using Gridap.ReferenceFEs: _CDLagrangianRefFE

T = Float64
reffe = _CDLagrangianRefFE(T,SEGMENT,(2,),(DISC,))
test_lagrangian_reference_fe(reffe)

T = Float64
reffe = _CDLagrangianRefFE(T,QUAD,1,(CONT,DISC))
test_lagrangian_reference_fe(reffe)

conf = CDConformity( (CONT,DISC) )
@test conf == Conformity(reffe)
@test get_face_own_dofs(QUAD4,conf) == get_face_own_dofs(reffe)
@test get_face_own_dofs_permutations(QUAD4,conf) == get_face_own_dofs_permutations(reffe)

T = VectorValue{2,Float64}
reffe = _CDLagrangianRefFE(T,QUAD,(2,0),(CONT,DISC))
test_lagrangian_reference_fe(reffe)

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4, 5, 6]]
@test get_face_own_nodes(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3]]

T = VectorValue{3,Float64}
reffe = _CDLagrangianRefFE(T,HEX,(2,2,2),(CONT,CONT,DISC))
test_lagrangian_reference_fe(reffe)

T = VectorValue{3,Float64}
reffe = _CDLagrangianRefFE(T,HEX,2,(CONT,CONT,DISC))
test_lagrangian_reference_fe(reffe)

T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,QUAD,(2,0))
@test Conformity(reffe) == CDConformity((CONT,DISC))
test_lagrangian_reference_fe(reffe)

T = VectorValue{3,Float64}
reffe = LagrangianRefFE(T,HEX,(0,2,0))
@test Conformity(reffe) == CDConformity((DISC,CONT,DISC))
test_lagrangian_reference_fe(reffe)

reffe = ReferenceFE(HEX,lagrangian,T,(0,2,0))
@test Conformity(reffe) == CDConformity((DISC,CONT,DISC))
test_lagrangian_reference_fe(reffe)

reffe = LagrangianRefFE(T,QUAD,(2,2))

end # module
