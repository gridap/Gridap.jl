# module CDLagrangianRefFEsTests

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs

T = Float64
reffe = CDLagrangianRefFE(T,SEGMENT,(2,),(DISC,))
test_reference_fe(reffe)

T = VectorValue{2,Float64}
reffe = CDLagrangianRefFE(T,QUAD,(2,0),(CONT,DISC))
test_reference_fe(reffe)

T = VectorValue{3,Float64}
reffe = CDLagrangianRefFE(T,HEX,(2,2,2),(CONT,CONT,DISC))
test_reference_fe(reffe)

# end # module
