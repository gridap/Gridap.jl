module PDiscRefFEsTests

using Gridap.TensorValues
using Gridap.ReferenceFEs

T = VectorValue{2,Float64}
reffe = PDiscRefFE(T,QUAD,2)
test_nodal_reference_fe(reffe)

T = VectorValue{2,Float64}
reffe = PDiscRefFE(T,HEX,3)
test_nodal_reference_fe(reffe)

end # module
