module MappingsTests

using Test
using Gridap.Arrays: CachedArray
using Gridap.Mappings
using Gridap.TensorValues
# using LinearAlgebra
# using Gridap.Inference

test_mapping(FunctionMapping(+),(3,2),5)
@test Mappings.return_types((FunctionMapping(+),FunctionMapping(/)),1,1) == (Int,Float64)

test_mapping(+,(3,2),5)

f = BroadcastMapping(+)
a = rand(3,2)
b = 3
c = a .+ b
Mappings.test_mapping(f,(a,b),c)


k = BroadcastMapping(-)
test_mapping(k,(1,),-1)
test_mapping(k,([1,2],),[-1,-2])
test_mapping(k,(1,2),-1)
test_mapping(k,(1.0,2),-1.0)
test_mapping(k,(1,2.0),-1.0)
test_mapping(k,([1,2],2),[-1,0])
test_mapping(k,(2,[1,2]),[1,0])
test_mapping(k,([3,4],[1,2]),[2,2])

f = BroadcastMapping(⋅)
a = fill(TensorValue(2,0,0,0,2,0,0,0,2),2)
b = VectorValue(1,2,3)
c = zeros(VectorValue{3,Int},2)
broadcast!(⋅,c,a,b)
test_mapping(f,(a,b),c)

end # module
