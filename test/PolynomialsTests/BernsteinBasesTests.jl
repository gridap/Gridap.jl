module BernsteinBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.Polynomials: binoms
using ForwardDiff

@test isHierarchical(Bernstein) == false

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
b = BernsteinBasis(Val(1),V,order)
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
b = BernsteinBasis(Val(1),V,order)

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
b = BernsteinBasis(Val(1),V,order)

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

function bernstein(K,N)
  b = binoms(Val(K))
  t -> b[N+1]*(t^N)*((1-t)^(K-N))
end
_∇(b) = t -> ForwardDiff.derivative(b, t)
_H(b) = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(b, y), t)

order = 3
b = BernsteinBasis(Val(1),V,order)

# x=x^1; x2 = x^2; x3 = x^3
#    -x3+3x2-3x+1  3x3-6x2+3x -3x3+3x2 x3
bx  = [      bernstein(order,n)( xi[1])  for xi in x,  n in 0:order]

#   -3x2+6x-3  9x2-12x+3  -9x2+6x 3x2
∇bx = [ G(_∇(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]

#    -6x+6 18x-12 -18x+6  6x
Hbx = [ H(_H(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


# Order 4

order = 4
b = BernsteinBasis(Val(1),V,order)

bx  = [      bernstein(order,n)( xi[1])  for xi in x,  n in 0:order]
∇bx = [ G(_∇(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]
Hbx = [ H(_H(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


end # module
