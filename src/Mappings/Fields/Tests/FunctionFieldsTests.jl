module FunctionFieldsTests

using Gridap.Inference
using Gridap.Mappings
using Gridap.NewFields
using Gridap.Arrays: CachedArray

using Gridap.TensorValues

u(x) = x


p1 = Point(1,2)
p2 = Point(2,1)
p3 = Point(4,3)
p4 = Point(6,1)
x = [p1,p2,p3,p4]

f = FunctionField(u)

test_field(f,x,CachedArray(x))

test_mapping(u,(x,),x)

end #module
