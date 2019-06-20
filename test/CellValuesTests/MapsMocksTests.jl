module MapsMocksTests

using Test
using Gridap
using TensorValues

using ..MapsMocks

a = VectorValue(10,10)
b = VectorValue(15,20)
p1 = VectorValue(1,1)
p2 = VectorValue(2,2)
p3 = VectorValue(3,3)
p = [p1,p2,p3]
ao = [  a+pj for pj in p  ]
m = MockMap(a)
test_map(m,p,ao)

ao = [p[i]*j+a  for j in 1:3, i in 1:length(p)]
m = TestMap(a,3)
test_map(m,p,ao)

a = VectorValue(10.1,10.3)
b = VectorValue(15.3,20.5)
ao = [  a+pj for pj in p  ]
m = MockMap(a)
test_map(m,p,ao)

ao = [p[i]*j+a  for j in 1:3, i in 1:length(p)]
m = TestMap(a,3)
test_map(m,p,ao)

end # module
