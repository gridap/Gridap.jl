module CDLagrangianRefFEsTests

using Gridap.TensorValues
using Gridap.ReferenceFEs

T = VectorValue{2,Float64}
reffe = CDLagrangianRefFE(T,QUAD,2,(CONT,DISC))
test_nodal_reference_fe(reffe)

T = VectorValue{2,Float64}
reffe = CDLagrangianRefFE(T,HEX,3,(CONT,CONT,DISC))
test_nodal_reference_fe(reffe)

end # module
