module BernsteinBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Polynomials: binoms
using ForwardDiff

@test isHierarchical(Bernstein) == false

x = [Point(0.),Point(1.),Point(.4)]
x1 = x[1]

#####################################
# Tests for 1D Bernstein polynomial #
#####################################

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

function test_internals(order,x,bx,∇bx,Hbx)
  sz = (1,order+1)
  for (i,xi) in enumerate(x)
    v2 = zeros(sz)
    Polynomials._evaluate_1d!(Bernstein,Val(order),v2,xi,1)
    @test all( [ bxi[1]≈vxi[1] for (bxi,vxi) in zip(bx[i,:],v2[:,1]) ] )

    g2 = zeros(sz)
    Polynomials._gradient_1d!(Bernstein,Val(order),g2,xi,1)
    @test all( [ bxi[1]≈vxi[1] for (bxi,vxi) in zip(∇bx[i,:],g2[:,1]) ] )

    h2 = zeros(sz)
    Polynomials._hessian_1d!(Bernstein,Val(order),h2,xi,1)
    @test all( [ bxi[1]≈vxi[1] for (bxi,vxi) in zip(Hbx[i,:],h2[:,1]) ] )
  end
end


# order 0 degenerated case

order = 0
b = BernsteinBasis(Val(1),V,order)
@test get_order(b) == 0
@test get_orders(b) == (0,)

bx =   [ 1.; 1.; 1.;; ]
∇bx = G[ 0.; 0.; 0.;; ]
Hbx = H[ 0.; 0.; 0.;; ]

test_internals(order,x,bx,∇bx,Hbx)

test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


# Order 1

order = 1
b = BernsteinBasis(Val(1),V,order)

bx = [ 1.0  0.
       0.  1.0
       0.6  0.4]

∇bx = G[ -1. 1.
         -1. 1.
         -1. 1. ]

Hbx = H[  0. 0.
          0. 0.
          0. 0. ]

test_internals(order,x,bx,∇bx,Hbx)

test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


# Order 2

order = 2
b = BernsteinBasis(Val(1),V,order)

bx  = [ 1.0   0.    0.
        0.    0.    1.0
        0.36  0.48  0.16]


∇bx = G[ -2.  2.   0.
          0. -2.   2.
        -1.2  0.4  .8 ]

Hbx = H[  2. -4. 2.
          2. -4. 2.
          2. -4. 2. ]

test_internals(order,x,bx,∇bx,Hbx)

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


# Order 3

function bernstein(K,N)
  b = binoms(Val(K))
  t -> b[N+1]*(t^N)*((1-t)^(K-N))
end
_∇(b) = t -> ForwardDiff.derivative(b, t)
_H(b) = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(b, y), t)

_bx_1D( order,x)   = [      bernstein(order,n)( xi[1])  for xi in x,  n in 0:order]
_∇bx_1D(order,x,G) = [ G(_∇(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]
_Hbx_1D(order,x,H) = [ H(_H(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]

order = 3
b = BernsteinBasis(Val(1),V,order)

# x=x^1; x2 = x^2; x3 = x^3
#    -x3+3x2-3x+1  3x3-6x2+3x -3x3+3x2 x3
bx  = _bx_1D( order,x)
#   -3x2+6x-3  9x2-12x+3  -9x2+6x 3x2
∇bx = _∇bx_1D(order,x,G)
#    -6x+6 18x-12 -18x+6  6x
Hbx = _Hbx_1D(order,x,H)
test_internals(order,x,bx,∇bx,Hbx)

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


# Order 4

order = 4
b = BernsteinBasis(Val(1),V,order)

bx  = _bx_1D( order,x)
∇bx = _∇bx_1D(order,x,G)
Hbx = _Hbx_1D(order,x,H)

test_internals(order,x,bx,∇bx,Hbx)

test_field_array(b,x,bx,≈, grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])


#####################################
# Tests for ND Bernstein polynomial #
#####################################

function bernstein_nD(D,K)
  indexbase = 1
  terms = Polynomials.bernstein_terms(Val(K),Val(D))
  coefs = Polynomials.multinoms(Val(K),Val(D))
  N = length(coefs)

  x -> begin
    @assert length(x) == D
    vals = zeros(eltype(x),N)

    # change to barycentric coords
    λ = zeros(eltype(x), D+1)
    copyto!(λ, x)
    λ[D+1] = 1 - sum(x)

    # apply formula
    for i in 1:N
      vals[i] = coefs[i]
      for (λi,ei) in zip(λ,terms[i])
        vals[i] *= λi^(ei-indexbase)
      end
    end
    return vals
  end
end

_∇(b) = x -> ForwardDiff.jacobian(b, get_array(x))
_H(b) = x -> ForwardDiff.jacobian(y -> ForwardDiff.jacobian(b, y), get_array(x))

_bx( D,order,x)   = transpose(reduce(hcat, (                                  bernstein_nD(D,order )(xi)           for xi in x)))
_∇bx(D,order,x,G) = transpose(reduce(hcat, ( map(G,        eachrow(        _∇(bernstein_nD(D,order))(xi)))         for xi in x)))
_Hbx(D,order,x,H) = transpose(reduce(hcat, ( map(splat(H), eachrow(reshape(_H(bernstein_nD(D,order))(xi), :,D*D))) for xi in x)))


D = 2
x = [Point(1.,0.), Point(.0,.5), Point(1.,.5), Point(.2,.3)]
x1 = x[1]

# scalar valued in 2D
V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 0
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
∇bx = _∇bx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

order = 1
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
∇bx = _∇bx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

order = 2
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
∇bx = _∇bx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

order = 3
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
∇bx = _∇bx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

order = 4
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
∇bx = _∇bx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)


# Vector valued in 2D
V = VectorValue{3,Float64}
G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 2
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = V[(1., 0., 0.) (0., 1., 0.) (0., 0., 1.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.);
        (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25) (0.5, 0., 0.) (0., 0.5, 0.) (0., 0., 0.5) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25);
        (1., 0., 0.) (0., 1., 0.) (0., 0., 1.) (1., 0., 0.) (0., 1., 0.) (0., 0., 1.) (-1., 0., 0.) (0., -1., 0.) (0., 0., -1.) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25) (-0.5, 0., 0.) (0., -0.5, 0.) (0., 0., -0.5) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25);
        (0.04, 0., 0.) (0., 0.04, 0.) (0., 0., 0.04) (0.12, 0., 0.) (0., 0.12, 0.) (0., 0., 0.12) (0.2, 0., 0.) (0., 0.2, 0.) (0., 0., 0.2) (0.09, 0., 0.) (0., 0.09, 0.) (0., 0., 0.09) (0.3, 0., 0.) (0., 0.3, 0.) (0., 0., 0.3) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25)]
∇bx = G[(2., 0., 0., 0., 0., 0.) (0., 0., 2., 0., 0., 0.) (0., 0., 0., 0., 2., 0.) (0., 2., 0., 0., 0., 0.) (0., 0., 0., 2., 0., 0.) (0., 0., 0., 0., 0., 2.) (-2., -2., 0., 0., 0., 0.) (0., 0., -2., -2., 0., 0.) (0., 0., 0., 0., -2., -2.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (-0., 0., 0., 0., 0., 0.) (0., 0., -0., 0., 0., 0.) (0., 0., 0., 0., -0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.);
        (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (1., 0., 0., 0., 0., 0.) (0., 0., 1., 0., 0., 0.) (0., 0., 0., 0., 1., 0.) (1., -0., 0., 0., 0., 0.) (0., 0., 1., -0., 0., 0.) (0., 0., 0., 0., 1., -0.) (0., 1., 0., 0., 0., 0.) (0., 0., 0., 1., 0., 0.) (0., 0., 0., 0., 0., 1.) (-1., 0., 0., 0., 0., 0.) (0., 0., -1., 0., 0., 0.) (0., 0., 0., 0., -1., 0.) (-1., -1., 0., 0., 0., 0.) (0., 0., -1., -1., 0., 0.) (0., 0., 0., 0., -1., -1.);
        (2., 0., 0., 0., 0., 0.) (0., 0., 2., 0., 0., 0.) (0., 0., 0., 0., 2., 0.) (1., 2., 0., 0., 0., 0.) (0., 0., 1., 2., 0., 0.) (0., 0., 0., 0., 1., 2.) (-3., -2., 0., 0., 0., 0.) (0., 0., -3., -2., 0., 0.) (0., 0., 0., 0., -3., -2.) (0., 1., 0., 0., 0., 0.) (0., 0., 0., 1., 0., 0.) (0., 0., 0., 0., 0., 1.) (-1., -2., 0., 0., 0., 0.) (0., 0., -1., -2., 0., 0.) (0., 0., 0., 0., -1., -2.) (1., 1., 0., 0., 0., 0.) (0., 0., 1., 1., 0., 0.) (0., 0., 0., 0., 1., 1.);
        (0.4, 0., 0., 0., 0., 0.) (0., 0., 0.4, 0., 0., 0.) (0., 0., 0., 0., 0.4, 0.) (0.6, 0.4, 0., 0., 0., 0.) (0., 0., 0.6, 0.4, 0., 0.) (0., 0., 0., 0., 0.6, 0.4) (0.6, -0.4, 0., 0., 0., 0.) (0., 0., 0.6, -0.4, 0., 0.) (0., 0., 0., 0., 0.6, -0.4) (0., 0.6, 0., 0., 0., 0.) (0., 0., 0., 0.6, 0., 0.) (0., 0., 0., 0., 0., 0.6) (-0.6, 0.4, 0., 0., 0., 0.) (0., 0., -0.6, 0.4, 0., 0.) (0., 0., 0., 0., -0.6, 0.4) (-1., -1., 0., 0., 0., 0.) (0., 0., -1., -1., 0., 0.) (0., 0., 0., 0., -1., -1.)]
Hbx = H[(2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (-4., -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4., -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4., -2., -2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.) (0., -2., -2., -4., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.) (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.);
        (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (-4., -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4., -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4., -2., -2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.) (0., -2., -2., -4., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.) (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.);
        (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (-4., -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4., -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4., -2., -2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.) (0., -2., -2., -4., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.) (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.);
        (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (-4., -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4., -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4., -2., -2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.) (0., -2., -2., -4., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.) (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.)]
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

# scalar valued in 3D

x = [Point(0.,0.,1.), Point(.5,.5,.5), Point(1.,.2,.4), Point(.2,.4,.3)]
x1 = x[1]

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 4
D = 3
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
∇bx = _∇bx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

b = BernsteinBasisOnSimplexDC(Val(D),V,order)
test_field_array(b,x,bx,≈, grad=∇bx, gradgrad=Hbx)

end # module
