module CurlConformBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials

PT = Monomial
T = Float64
k = 1

#################################################
# Curl conform bases on n-cubes (non Bernstein) #
#################################################

xi = Point(2,3)
np = 5
x = fill(xi,np)

# 3D
order = 0
r = order+1
D = 2
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT)

@test length(b) == 4
@test get_order(b) == 1
@test testvalue(typeof(b)) isa typeof(b)

V = VectorValue{D,T}
@test_throws AssertionError FEEC_poly_basis(Val(D),V,r,k,:Q⁻,PT)

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 0
r = order+1
D = 3
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
H = gradient_type(G,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT)

v = V[
# (1.0, 0.0, 0.0), ( y , 0.0, 0.0), ( z , 0.0, 0.0), ( yz , 0.0, 0.0),
# (0.0, 1.0, 0.0), (0.0,  x , 0.0), (0.0,  z , 0.0), (0.0,  xz , 0.0),
# (0.0, 0.0, 1.0), (0.0, 0.0,  x ), (0.0, 0.0,  y ), (0.0, 0.0,  xy)]
  (1.0, 0.0, 0.0), (3.0, 0.0, 0.0), (5.0, 0.0, 0.0), (15.0, 0.0, 0.0),
  (0.0, 1.0, 0.0), (0.0, 2.0, 0.0), (0.0, 5.0, 0.0), (0.0, 10.0, 0.0),
  (0.0, 0.0, 1.0), (0.0, 0.0, 2.0), (0.0, 0.0, 3.0), (0.0, 0.0, 6.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 1 ,  0 , 0  )
  (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( y ,  0 , 0  )
  (0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( z ,  0 , 0  )
  (0.0, 5.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( yz , 0 , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  1 , 0  )
  (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  x , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0), # ( 0 ,  z , 0  )
  (0.0, 0.0, 0.0, 5.0, 0.0, 2.0, 0.0, 0.0, 0.0), # ( 0 ,  xz, 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  0 , 1  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0), # ( 0 ,  0 , x  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0), # ( 0 ,  0 , y  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 2.0, 0.0)] # ( 0 ,  0 , xy )

h = H[
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 1 ,  0 , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( y ,  0 , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( z ,  0 , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 1. , 0.0, 1. , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( yz , 0 , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  1 , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  x , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  z , 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1. , 0.0, 0.0, 0.0, 1. , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  xz, 0  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  0 , 1  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  0 , x  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), # ( 0 ,  0 , y  )
  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1. , 0.0, 1. , 0.0, 0.0, 0.0, 0.0, 0.0)] # ( 0 ,  0 , xy )

bx  = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)

evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# 2D

xi = Point(2,3)
np = 5
x = fill(xi,np)

order = 1
r = order+1
D = 2
V = VectorValue{D,T}
G = gradient_type(V,xi)
H = gradient_type(G,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT)
@test testvalue(typeof(b)) isa typeof(b)

v = V[
# (1., 0.), (x , 0.), (y , 0.), (xy, 0.), (y², 0.), (xy², 0.),
# (0., 1.), (0., x ), (0., x²), (0., y ), (0., xy), (0., x²y)]
  (1., 0.), (2., 0.), (3., 0.), (6., 0.), (9., 0.), (18., 0.),
  (0., 1.), (0., 2.), (0., 4.), (0., 3.), (0., 6.), (0., 12.)]

g = G[
  #                    # ∂xV₁ ∂yV₁ ∂ₓV₂ ∂xV₂  # (V₁, V₂ )
  (0.,  0.,  0.,  0.), # (0 ,  0 ,  0 , 0 )  # (1 , 0  )
  (1.,  0.,  0.,  0.), # (1 ,  0 ,  0 , 0 )  # (x , 0  )
  (0.,  1.,  0.,  0.), # (0 ,  1 ,  0 , 0 )  # (y , 0  )
  (3.,  2.,  0.,  0.), # (y ,   x,  0 , 0 )  # (xy, 0  )
  (0.,  6.,  0.,  0.), # (0.,  2y,  0 , 0 )  # (y², 0  )
  (9., 12.,  0.,  0.), # (y², 2xy,  0 , 0 )  # (xy²,0  )
  (0.,  0.,  0.,  0.), # (0 ,  0 ,  0 , 0 )  # (0 , 1  )
  (0.,  0.,  1.,  0.), # (0 ,  0 ,  1 , 0 )  # (0 , x  )
  (0.,  0.,  4.,  0.), # (0 ,  0 ,  2x, 0 )  # (0 , x² )
  (0.,  0.,  0.,  1.), # (0 ,  0 ,  0 , 1 )  # (0 , y  )
  (0.,  0.,  3.,  2.), # (0 ,  0 ,   y, x )  # (0 , xy )
  (0.,  0., 12.,  4.)] # (0 ,  0 , 2xy, x²)  # (0 , x²y)

h = H[
 # ∂xxV₁ ∂yxV₁ ∂xyV₁ ∂yyV₁ ∂xxV₂ ∂yxV₂ ∂xxV₂ ∂yxV₂  # (V₁, V₂ )
  (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 ),  # (1., 0. ) # (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 ),  # (x , 0. ) # (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 ),  # (y , 0. ) # (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 )
  (0 ,  1 ,  1 ,  0 ,  0 ,  0 , 0 , 0 ),  # (xy, 0. ) # (0 ,  1 ,   1,   0,  0 ,  0 , 0 , 0 )
  (0 ,  0.,  0 ,  2 ,  0 ,  0 , 0 , 0 ),  # (y², 0. ) # (0 ,  0.,   0,   2,  0 ,  0 , 0 , 0 )
  (0 ,  6.,  6.,  4.,  0 ,  0 , 0 , 0 ),  # (xy², 0.) # (0 ,  2y,  2y,  2x,  0 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 ),  # (0., 1. ) # (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 ),  # (0., x  ) # (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  2 ,  0 , 0 , 0 ),  # (0., x² ) # (0 ,  0 ,  0 ,  0 ,  2 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 ),  # (0., y  ) # (0 ,  0 ,  0 ,  0 ,  0 ,  0 , 0 , 0 )
  (0 ,  0 ,  0 ,  0 ,  0 ,  1 , 1 , 0 ),  # (0., xy ) # (0 ,  0 ,  0 ,  0 ,  0 ,  1 , 1 , 0 )
  (0 ,  0 ,  0 ,  0 ,  6.,  4., 4., 0 )]  # (0., x²y) # (0 ,  0 ,  0 ,  0 ,  2y,  2x, 2x, 0 )

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)
test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


# 1D

order = 0
r = order+1
D = 1
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT)

@test b isa CartProdPolyBasis{D,T,PT}
@test testvalue(typeof(b)) isa typeof(b)

# Misc

# Derivatives not implemented for symetric tensor types

V = SymTensorValue{D,T}
G = gradient_type(V,xi)
r = zeros(G, (1,1))
@test_throws ErrorException Polynomials._comp_wize_set_derivative!(r,0,0,0,V)

V = SymTracelessTensorValue{D,T}
G = gradient_type(V,xi)
r = zeros(G, (1,1))
@test_throws ErrorException Polynomials._comp_wize_set_derivative!(r,0,0,0,V)


###################################################
# Curl conform bases on simplices (non Bernstein) #
###################################################

xi = Point(2.,3.,5.)
np = 3
x = fill(xi,np)

order = 0
D = 3
b = NedelecPolyBasisOnSimplex{D}(PT, T, order)
@test testvalue(typeof(b)) isa typeof(b)

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
@test testvalue(typeof(b)) isa typeof(b)
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
r = order+1
D = 1
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,PT)

@test b isa CartProdPolyBasis{D,T,PT}
@test testvalue(typeof(b)) isa typeof(b)

# D > 3 not implemented

@test_throws ErrorException NedelecPolyBasisOnSimplex{4}(PT, T, order)

end # module
