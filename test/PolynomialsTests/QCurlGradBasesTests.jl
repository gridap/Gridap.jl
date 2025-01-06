module QCurlGradBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

xi = Point(2,3)
np = 5
x = fill(xi,np)

PT = Monomial

order = 0
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = QCurlGradBasis(PT,Val(D),T,order)

@test length(b) == 4
@test get_order(b) == 1

@test_throws AssertionError QCurlGradBasis(PT,Val(D),V,order)

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 0
D = 3
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = QCurlGradBasis(PT, Val(D),T,order)

v = V[
  (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
  (0.0, 1.0, 0.0), (0.0, 3.0, 0.0),
  (0.0, 0.0, 1.0), (0.0, 0.0, 5.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
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
b = QCurlGradBasis(PT, Val(D),T,order)

v = V[
  ( 1., 0. ), ( 2., 0. ), ( 4., 0. ), ( 3., 0. ), ( 6., 0. ), (12., 0. ),
  ( 0., 1. ), ( 0., 2. ), ( 0., 3. ), ( 0., 6. ), ( 0., 9. ), ( 0., 18.)]

g = G[
  ( 0., 0., 0., 0.),
  ( 1., 0., 0., 0.),
  ( 4., 0., 0., 0.),
  ( 0., 1., 0., 0.),
  ( 3., 2., 0., 0.),
  (12., 4., 0., 0.),
  ( 0., 0., 0., 0.),
  ( 0., 0., 1., 0.),
  ( 0., 0., 0., 1.),
  ( 0., 0., 3., 2.),
  ( 0., 0., 0., 6.),
  ( 0., 0., 9.,12.)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

# 1D

order = 0
D = 1
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = QCurlGradBasis(PT,Val(D),T,order)

@test b isa UniformPolyBasis{D,V,order+1,PT}

@test_throws AssertionError QCurlGradBasis(PT,Val(D),V,order)


end # module
