# module CDLagrangianRefFEsTests

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Test

T = Float64
reffe = CDLagrangianRefFE(T,SEGMENT,(2,),(DISC,))
test_reference_fe(reffe)

T = Float64
reffe = CDLagrangianRefFE(T,QUAD,1,(CONT,DISC))
test_reference_fe(reffe)

conf = CDConformity( (CONT,DISC) )
@test conf == get_default_conformity(reffe)
@test get_face_own_dofs(QUAD4,conf) == get_face_own_dofs(reffe)
@test get_face_own_dofs_permutations(QUAD4,conf) == get_face_own_dofs_permutations(reffe)

T = VectorValue{2,Float64}
reffe = CDLagrangianRefFE(T,QUAD,(2,0),(CONT,DISC))
test_reference_fe(reffe)

T = VectorValue{3,Float64}
reffe = CDLagrangianRefFE(T,HEX,(2,2,2),(CONT,CONT,DISC))
test_reference_fe(reffe)

T = VectorValue{3,Float64}
reffe = CDLagrangianRefFE(T,HEX,2,(CONT,CONT,DISC))
test_reference_fe(reffe)

T = VectorValue{2,Float64}
reffe = LagrangianRefFE(T,QUAD,(2,0))
test_reference_fe(reffe)

reffe = LagrangianRefFE(T,QUAD,(2,2))

# end # module
