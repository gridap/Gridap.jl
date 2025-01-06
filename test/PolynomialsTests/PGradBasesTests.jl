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

PT = Monomial

order = 0
D = 3
b = NedelecPolyBasisOnSimplex{D}(PT, T, order)

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
b = NedelecPolyBasisOnSimplex{D}(PT, T, order)
V = VectorValue{D, Float64}
v = V[(1,0),(0,1),(-3,2)]
G = gradient_type(V,xi)
g = G[(0,0, 0,0), (0,0, 0,0), (0,-1, 1,0)]
bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

# 1D

order = 0
D = 1
T = Float64
V = VectorValue{D,T}
b = QGradBasis(PT,Val(D),T,order)

@test b isa UniformPolyBasis{D,V,order+1,PT}

@test_throws AssertionError QGradBasis(PT,Val(D),V,order)

# D > 3 not implemented

@test_throws ErrorException NedelecPolyBasisOnSimplex{4}(PT, T, order)

end # module
