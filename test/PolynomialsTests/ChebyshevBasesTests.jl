module ChebyshevBasisTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using ForwardDiff

@test isHierarchical(Chebyshev) == true

np = 3
x = [Point(2.),Point(3.),Point(4.)]
x1 = x[1]


# Only test 1D evaluations as tensor product structure is tested in monomial tests

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

chebyshev_T(N) = t -> begin
  ξ = 2t-1
  sq = sqrt(ξ*ξ-1)
  .5*( (ξ - sq)^N + (ξ + sq)^N )
end
chebyshev_U(N) = t -> begin
  ξ = 2t-1
  sq = sqrt(ξ*ξ-1)
  .5*( (ξ + sq)^(N+1) - (ξ - sq)^(N+1) )/sq
end
_∇(b) = t -> ForwardDiff.derivative(b, t)
_H(b) = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(b, y), t)

######################################
# First Kind Chebyshev ( Hessian TODO)
######################################

# order 0 degenerated case
order = 0
b = ChebyshevBasis(Val(1),V,order)
@test get_order(b) == 0
@test get_orders(b) == (0,)
@test testvalue(typeof(b)) isa typeof(b)

bx  = [      chebyshev_T(n)( xi[1])  for xi in x,  n in 0:order]
∇bx = [ G(_∇(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]
Hbx = [ H(_H(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈, grad=∇bx)#,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈,grad=∇bx[1,:])#,gradgrad=Hbx[1,:])


# Order 1
order = 1
b = ChebyshevBasis(Val(1),V,order)

bx  = [      chebyshev_T(n)( xi[1])  for xi in x,  n in 0:order]
∇bx = [ G(_∇(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]
Hbx = [ H(_H(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈,grad=∇bx)#,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈,grad=∇bx[1,:])#,gradgrad=Hbx[1,:])


# Order 2
order = 2
b = ChebyshevBasis(Val(1),V,order)


bx  = [      chebyshev_T(n)( xi[1])  for xi in x,  n in 0:order]
∇bx = [ G(_∇(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]
Hbx = [ H(_H(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈,grad=∇bx)#,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈,grad=∇bx[1,:])#,gradgrad=Hbx[1,:])


# Order 3
order = 3
b = ChebyshevBasis(Val(1),V,order)

bx  = [      chebyshev_T(n)( xi[1])  for xi in x,  n in 0:order]
∇bx = [ G(_∇(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]
Hbx = [ H(_H(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈, grad=∇bx)#,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈,grad=∇bx[1,:])#,gradgrad=Hbx[1,:])


# Order 4
order = 4
b = ChebyshevBasis(Val(1),V,order)

bx  = [      chebyshev_T(n)( xi[1])  for xi in x,  n in 0:order]
∇bx = [ G(_∇(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]
Hbx = [ H(_H(chebyshev_T(n))(xi[1])) for xi in x,  n in 0:order]

test_field_array(b,x,bx,≈,grad=∇bx)#,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈,grad=∇bx[1,:])#,gradgrad=Hbx[1,:])


############################
# Second Kind Chebyshev TODO
############################

@test_throws ErrorException ChebyshevBasis(Val(1),Float64,0;kind=:U)

order = 1
d = 1
b = ChebyshevBasis(Val(1),V,order)
Hb = FieldGradientArray{2}(b)
r, c, g, h = return_cache(Hb,x)

@test_throws ErrorException Polynomials._gradient_1d!(Chebyshev{:U},Val(order),g,x1,d)
@test_throws ErrorException Polynomials._hessian_1d!( Chebyshev{:U},Val(order),h,x1,d)


end # module
