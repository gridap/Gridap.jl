module QCurlGradMonomialBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

xi = Point(2,3)
np = 5
x = fill(xi,np)

order = 0
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = QCurlGradMonomialBasis{D}(T,order)

@test num_terms(b) == 4
@test get_order(b) == 0

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 0
D = 3
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = QCurlGradMonomialBasis{D}(T,order)

v = V[
  (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
  (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 5.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

xi = Point(2,3)
np = 5
x = fill(xi,np)

order = 1
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = QCurlGradMonomialBasis{D}(T,order)

v = V[
  (1.0, 0.0), (0.0, 1.0), (2.0, 0.0), (0.0, 3.0),
  (4.0, 0.0), (0.0, 9.0), (3.0, 0.0), (0.0, 2.0),
  (6.0, 0.0), (0.0, 6.0), (12.0, 0.0), (0.0, 18.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0),
  (1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0),
  (4.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 6.0),
  (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0),
  (3.0, 2.0, 0.0, 0.0), (0.0, 0.0, 3.0, 2.0),
  (12.0, 4.0, 0.0, 0.0),(0.0, 0.0, 9.0, 12.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

end # module
