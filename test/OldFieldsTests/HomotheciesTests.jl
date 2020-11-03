module HomotheciesTests

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields

origin = Point(1,1)
scaling = Point(2,2)

h = Homothecy(origin,scaling)

x1 = Point(0,0)
x2 = Point(1,1)
x3 = Point(2,1)
x = [x1,x2,x3]
r = Point{2,Int}[(-1, -1), (1, 1), (3, 1)]
∇r = TensorValue{2,2,Int,4}[(2, 0, 0, 2), (2, 0, 0, 2), (2, 0, 0, 2)]
test_field(h,x,r,grad=∇r)

end # module
