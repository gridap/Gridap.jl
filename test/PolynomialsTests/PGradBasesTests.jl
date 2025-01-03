module PGradBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials

xi = Point(2.,3.,5.)
np = 3
x = fill(xi,np)
T = Float64


order = 0
D = 3
b = NedelecPolyBasisOnSimplex{D}(Monomial, T, order)

V = VectorValue{D, Float64}
v = V[(1,0,0),(0,1,0),(0,0,1),(-3,2,0),(-5,0,2),(0,-5,3)]

G = gradient_type(V,xi)
g = G[
  (0,0,0, 0,0,0, 0,0,0), (0,0,0, 0,0,0, 0,0,0), (0,0,0, 0,0,0, 0,0,0),
  (0,-1,0, 1,0,0, 0,0,0),(0,0,-1, 0,0,0, 1,0,0),(0,0,0, 0,0,-1, 0,1,0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

xi = Point(2.,3.)
np = 4
x = fill(xi,np)

order = 0
D = 2
b = NedelecPolyBasisOnSimplex{D}(Monomial, T, order)
V = VectorValue{D, Float64}
v = V[(1,0),(0,1),(-3,2)]
G = gradient_type(V,xi)
g = G[(0,0, 0,0), (0,0, 0,0), (0,-1, 1,0)]
bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])


end # module
