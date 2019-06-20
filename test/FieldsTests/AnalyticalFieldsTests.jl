module AnalyticalFieldsTests

using Test
using Gridap
import Gridap: ∇

p1 = Point(1.0,2.0)
p2 = Point(0.0,3.0)
p3 = Point(7.0,4.0)
p4 = Point(8.0,1.0)
p = [p1,p2,p3,p4]

fun(x) = 3*x[1]

f = AnalyticalField(fun,2)
v = fun.(p)
test_field_without_gradient(f,p,v)

fungrad(x) = VectorValue(3.0,0.0)

∇(::typeof(fun)) = fungrad

f = AnalyticalField(fun,2)
v = fun.(p)
g = fungrad.(p)
test_field(f,p,v,g)

end # module
