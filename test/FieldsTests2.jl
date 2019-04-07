module FieldsTests2

using Test
using Numa.Fields2
using Numa.FieldValues

import Numa.Fields2: gradient

fun(x::Point{2}) = x[1]*x[2] + x[1]

gradfun(x::Point{2}) = VectorValue(x[2]+1.0,x[1])

gradgradfun(x::Point{2}) = TensorValue(0.0,1.0,1.0,0.0)

gradient(::typeof(fun)) = gradfun

gradient(::typeof(gradfun)) = gradgradfun

p1 = VectorValue(0.1,1.0)
p2 = VectorValue(1.1,2.0)
p3 = VectorValue(1.4,5.0)

p = [p1,p2,p3]

f = AnalyticalField(fun,2)
@test evaluate(f,p) == fun.(p)

∇f = gradient(f)
@test evaluate(∇f,p) == gradfun.(p)

∇∇f = gradient(∇f)
@test evaluate(∇∇f,p) == gradgradfun.(p)

end #module FieldsTests2

