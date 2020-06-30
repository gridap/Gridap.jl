module CDLagrangianRefFEsTests

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs

T = Float64
reffe = CDLagrangianRefFE(T,SEGMENT,(2,),(DISC,))
test_reference_fe(reffe)

reffe = CDLagrangianRefFE(T,QUAD,(2,2),(CONT,DISC))
test_reference_fe(reffe)

T = VectorValue{3,Float64}
reffe = CDLagrangianRefFE(T,HEX,(2,2,2),(CONT,CONT,DISC))
test_reference_fe(reffe)

T = Float64
reffe = LagrangianRefFE(T,QUAD,2)
reffe.face_own_nodes

end # module
