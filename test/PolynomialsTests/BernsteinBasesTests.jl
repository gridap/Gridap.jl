module BernsteinBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

np = 3
x = [Point(0.),Point(1.),Point(.4)]
xi = x[1]

# Only test 1D evaluations as tensor product structure is tested in monomial tests
#
V = Float64
G = gradient_type(V,xi)
H = gradient_type(G,xi)

# order 0 degenerated case

order = 0
b = BernsteinBasis{1}(V,order)
@test get_order(b) == 0
@test get_orders(b) == (0,)

v = V[1.0,]
g = G[(0.0)]
h = H[(0.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)
test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Order 1

order = 1
b = BernsteinBasis{1}(V,order)

bx = [ 1.0  0.0
       0.0  1.0
       0.6  0.4]

∇bx = G[ -1. 1.
         -1. 1.
         -1. 1. ]

Hbx = H[  0. 0.
          0. 0.
          0. 0. ]

test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Order 2

order = 2
b = BernsteinBasis{1}(V,order)

bx  = [ 1.0   0.0   0.0
        0.0   0.0   1.0
        0.36  0.48  0.16]


∇bx = G[ -2.  2.   0.
          0. -2.   2.
        -1.2  0.4  .8 ]

Hbx = H[  2. -4. 2.
          2. -4. 2.
          2. -4. 2. ]

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Order 3

order = 3
b = BernsteinBasis{1}(V,order)

# x=x^1; x2 = x^2; x3 = x^3
#    -x3+3x2-3x+1  3x3-6x2+3x -3x3+3x2 x3
bx  = [ 1.0   0.0  0.0   0.0
        0.0   0.0  0.0   1.0
        .216  .432 .288  .064]

#   -3x2+6x-3  9x2-12x+3  -9x2+6x 3x2
∇bx = G[ -3.   3.   0. 0.
          0.   0.  -3. 3.
        -1.08 -.36 .96 .48]

#    -6x+6 18x-12 -18x+6  6x
Hbx = H[ 6.  -12.   6.   0.
         0.   6.   -12.  6.
         3.6 -4.8  -1.2  2.4]

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

end # module
