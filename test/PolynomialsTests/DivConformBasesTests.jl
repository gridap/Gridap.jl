module DivConformBasesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

PT = Monomial
T = Float64

################################################
# Div conform bases on n-cubes (non Bernstein) #
################################################

xi = Point(2,3)
np = 5
x = fill(xi,np)

order = 0
D = 2
k, r, rotate_90 = D-1, order+1, D==2
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT; rotate_90)

@test length(b) == 4
@test get_order(b) == 1
@test get_orders(b) == (1,1)
@test value_type(b) == V
@test testvalue(typeof(b)) isa typeof(b)

@test_throws AssertionError FEEC_poly_basis(Val(D),V,r,k,:Q⁻,PT; rotate_90)

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 0
D = 3
k, r, rotate_90 = D-1, order+1, D==2
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT; rotate_90)

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
k, r, rotate_90 = D-1, order+1, D==2
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT; rotate_90)

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
k, r, rotate_90 = D-1, order+1, D==2
b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,PT; rotate_90)

@test b isa CartProdPolyBasis{D,T,PT}
@test testvalue(typeof(b)) isa typeof(b)


##################################################
# Div conform bases on simplices (non Bernstein) #
##################################################

xi = Point(4,2)
np = 1
x = fill(xi,np)

order = 2
D = 2
k, r, rotate_90 = D-1, order+1, D==2
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,PT; rotate_90)

@test testvalue(typeof(b)) isa typeof(b)

v = V[
  (1.0,  0.0), (4.0,  0.0), (16.0, 0.0), (2.0,  0.0), (8.0,  0.0), (4.0,  0.0), # pterm ex
  (0.0,  1.0), (0.0,  4.0), (0.0, 16.0), (0.0,  2.0), (0.0,  8.0), (0.0,  4.0), # pterm ey
  (64.0,32.0), (32.0,16.0), (16.0, 8.0)] # sterms

g = G[
  (0.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0), (8.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (2.0, 4.0, 0.0, 0.0), (0.0, 4.0, 0.0, 0.0), # pterm ex
  (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (0.0, 0.0, 8.0, 0.0), (0.0, 0.0, 0.0, 1.0), (0.0, 0.0, 2.0, 4.0), (0.0, 0.0, 0.0, 4.0), # pterm ey
  (48.0, 0.0, 16.0, 16.0), (16.0, 16.0, 4.0, 16.0), (4.0, 16.0, 0.0, 12.0)] # sterms

vb = evaluate(b,x)

for (vi,vbi) in zip(v,vb)
  @test vi == vbi
end

vb = evaluate(b,xi)
@test vb == v

∇b = Broadcasting(gradient)(b)
gvb = evaluate(∇b,x)
for (vi,vbi) in zip(g,gvb)
  @test vi == vbi
end

gvb = evaluate(∇b,xi)
@test gvb == g

@test length(b) == 15
@test get_order(b) == 3
@test get_orders(b) == (3,3)

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 1
D = 3
k, r, rotate_90 = D-1, order+1, D==2
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,PT; rotate_90)

@test length(b) == 15
@test get_order(b) == 2
@test get_orders(b) == (2,2,2)


# 1D

order = 0
D = 1
k, r, rotate_90 = D-1, order+1, D==2
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,PT; rotate_90)
@test testvalue(typeof(b)) isa typeof(b)

@test b isa CartProdPolyBasis{D,T,PT}

end # module
