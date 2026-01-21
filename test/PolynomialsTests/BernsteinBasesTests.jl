module BernsteinBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using ForwardDiff
using StaticArrays

@test isHierarchical(Bernstein) == false

x = [Point(0.),Point(1.),Point(.4)]
x1 = x[1]

#####################################
# Tests for 1D Bernstein polynomial #
#####################################

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

function test_internals(order,x,bx,Gbx,Hbx)
  sz = (1,order+1)
  for (i,xi) in enumerate(x)
    v2 = zeros(sz)
    Polynomials._evaluate_1d!(Bernstein,order,v2,xi,1)
    @test all( [ bxi[1]≈vxi[1] for (bxi,vxi) in zip(bx[i,:],v2[:,1]) ] )

    g2 = zeros(sz)
    Polynomials._gradient_1d!(Bernstein,order,g2,xi,1)
    @test all( [ bxi[1]≈vxi[1] for (bxi,vxi) in zip(Gbx[i,:],g2[:,1]) ] )

    h2 = zeros(sz)
    Polynomials._hessian_1d!(Bernstein,order,h2,xi,1)
    @test all( [ bxi[1]≈vxi[1] for (bxi,vxi) in zip(Hbx[i,:],h2[:,1]) ] )
  end
end


# order 0 degenerated case

order = 0
b = BernsteinBasis(Val(1),V,order)
@test get_order(b) == 0
@test get_exponents(b) == [(0,),]
@test testvalue(typeof(b)) isa typeof(b)

bx =   [ 1.; 1.; 1.;; ]
Gbx = G[ (0.,); (0.,); (0.,);; ]
Hbx = H[ (0.,); (0.,); (0.,);; ]

test_internals(order,x,bx,Gbx,Hbx)

test_field_array(b,x,bx,grad=Gbx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])


# Order 1

order = 1
b = BernsteinBasis(Val(1),V,order)

bx = [ 1.0  0.
       0.  1.0
       0.6  0.4]

Gbx = G[ (-1.,) (1.,)
         (-1.,) (1.,)
         (-1.,) (1.,) ]

Hbx = H[  (0.,) (0.,)
          (0.,) (0.,)
          (0.,) (0.,) ]

test_internals(order,x,bx,Gbx,Hbx)

test_field_array(b,x,bx,grad=Gbx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])


# Order 2

order = 2
b = BernsteinBasis(Val(1),V,order)

bx  = [ 1.0   0.    0.
        0.    0.    1.0
        0.36  0.48  0.16]


Gbx = G[( -2.,) ( 2. ,)  (0.,)
        (  0.,) (-2. ,)  (2.,)
        (-1.2,) ( 0.4,)  (.8,) ]

Hbx = H[  (2.,) (-4.,) (2.,)
          (2.,) (-4.,) (2.,)
          (2.,) (-4.,) (2.,) ]

test_internals(order,x,bx,Gbx,Hbx)

test_field_array(b,x,bx,≈, grad=Gbx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])


# Order 3

function bernstein(K,N)
  b = ntuple( i -> binomial(K,i-1), Val(K+1)) # all binomial(i,K) for 0≤i≤K
  t -> b[N+1]*(t^N)*((1-t)^(K-N))
end
_∇(b) = t -> ForwardDiff.derivative(b, t)
_H(b) = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(b, y), t)

_bx_1D( order,x)   = [      bernstein(order,n)( xi[1])  for xi in x,  n in 0:order]
_Gbx_1D(order,x,G) = [ G(_∇(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]
_Hbx_1D(order,x,H) = [ H(_H(bernstein(order,n))(xi[1])) for xi in x,  n in 0:order]

order = 3
b = BernsteinBasis(Val(1),V,order)

# x=x^1; x2 = x^2; x3 = x^3
#    -x3+3x2-3x+1  3x3-6x2+3x -3x3+3x2 x3
bx  = _bx_1D( order,x)
#   -3x2+6x-3  9x2-12x+3  -9x2+6x 3x2
Gbx = _Gbx_1D(order,x,G)
#    -6x+6 18x-12 -18x+6  6x
Hbx = _Hbx_1D(order,x,H)
test_internals(order,x,bx,Gbx,Hbx)

test_field_array(b,x,bx,≈, grad=Gbx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])


# Order 4

order = 4
b = BernsteinBasis(Val(1),V,order)

bx  = _bx_1D( order,x)
Gbx = _Gbx_1D(order,x,G)
Hbx = _Hbx_1D(order,x,H)

test_internals(order,x,bx,Gbx,Hbx)

test_field_array(b,x,bx,≈, grad=Gbx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])


#####################################
# Tests for ND Bernstein polynomial #
#####################################

function bernstein_nD(D,K,x2λ=nothing)
  terms = Polynomials.bernstein_terms(K,D)
  N = length(terms)
  x -> begin
    @assert length(x) == D
    vals = zeros(eltype(x),N)
    # change to barycentric coords of reference simplex
    if isnothing(x2λ)
      λ = SVector{D+1,eltype(x)}(1 - sum(x), x...)
    else
      λ = x2λ*SVector(1, x...)
    end
    for i in 1:N
      vals[i] = Polynomials.multinomial(terms[i]...)
      for (λi,ei) in zip(λ,terms[i])
        vals[i] *= λi^ei
      end
    end
    return vals
  end
end

_∇(b) = x -> ForwardDiff.jacobian(b, get_array(x))
_H(b) = x -> ForwardDiff.jacobian(y -> ForwardDiff.jacobian(b, y), get_array(x))

_bx( D,order,x,  x2λ=nothing) = transpose(reduce(hcat, (                                    bernstein_nD(D,order,x2λ )(xi)           for xi in x)))
_Gbx(D,order,x,G,x2λ=nothing) = transpose(reduce(hcat, ( map(G,          eachrow(        _∇(bernstein_nD(D,order,x2λ))(xi)))         for xi in x)))
_Hbx(D,order,x,H,x2λ=nothing) = transpose(reduce(hcat, ( map(x->H(x...), eachrow(reshape(_H(bernstein_nD(D,order,x2λ))(xi), :,D*D))) for xi in x)))

D = 2
x = [Point(1.,0.), Point(.0,.5), Point(1.,.5), Point(.2,.3)]
x3 = x[3]

# scalar valued in 2D
V = Float64
G = gradient_type(V,x3)
H = gradient_type(G,x3)

order = 0
b = BernsteinBasisOnSimplex(Val(D),V,order)
@test get_order(b) == 0
bx  = _bx( D,order,x)
Gbx = _Gbx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)

order = 1
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
Gbx = _Gbx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)

order = 2
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
Gbx = _Gbx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)

order = 3
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
Gbx = _Gbx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)

order = 4
b = BernsteinBasisOnSimplex(Val(D),V,order)
bx  = _bx( D,order,x)
Gbx = _Gbx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],≈,grad=Gbx[1,:],gradgrad=Hbx[1,:])

b_terms = bernstein_terms(order,D)
λ = Polynomials._cart_to_bary(x3, b.cart_to_bary_matrix)
for j in 1:length(b)
  α = b_terms[j]
  id_α = bernstein_term_id(α)
  @test id_α == j
  c = Float64[ Int(i==id_α) for i in 1:length(b) ] # Bernstein coefficients of Bα
  Polynomials._de_Casteljau_nD!(c, λ, Val(order), Val(D))
  Bα_x3 = c[1]
  @test Bα_x3 == bx[3,id_α]
end

# Vector valued in 2D
V = VectorValue{3,Float64}
G = gradient_type(V,x3)
H = gradient_type(G,x3)

order = 2
b = BernsteinBasisOnSimplex(Val(D),V,order)
@test testvalue(typeof(b)) isa typeof(b)
bx  = V[(0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (1., 0., 0.) (0., 1., 0.) (0., 0., 1.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.);
        (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0.5, 0., 0.) (0., 0.5, 0.) (0., 0., 0.5) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0., 0., 0.) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25);
        (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25) (-1., 0., 0.) (0., -1., 0.) (0., 0., -1.) (-0.5, 0., 0.) (0., -0.5, 0.) (0., 0., -0.5) (1., 0., 0.) (0., 1., 0.) (0., 0., 1.) (1., 0., 0.) (0., 1., 0.) (0., 0., 1.) (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25); (0.25, 0., 0.) (0., 0.25, 0.) (0., 0., 0.25) (0.2, 0., 0.) (0., 0.2, 0.) (0., 0., 0.2) (0.3, 0., 0.) (0., 0.3, 0.) (0., 0., 0.3) (0.04000000000000001, 0., 0.) (0., 0.04000000000000001, 0.) (0., 0., 0.04000000000000001) (0.12, 0., 0.) (0., 0.12, 0.) (0., 0., 0.12) (0.09, 0., 0.) (0., 0.09, 0.) (0., 0., 0.09)]
Gbx = G[(0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (-2., -2., 0., 0., 0., 0.) (0., 0., -2., -2., 0., 0.) (0., 0., 0., 0., -2., -2.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (2., 0., 0., 0., 0., 0.) (0., 0., 2., 0., 0., 0.) (0., 0., 0., 0., 2., 0.) (0., 2., 0., 0., 0., 0.) (0., 0., 0., 2., 0., 0.) (0., 0., 0., 0., 0., 2.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.);
        (-1., -1., 0., 0., 0., 0.) (0., 0., -1., -1., 0., 0.) (0., 0., 0., 0., -1., -1.) (1., 0., 0., 0., 0., 0.) (0., 0., 1., 0., 0., 0.) (0., 0., 0., 0., 1., 0.) (-1., 0., 0., 0., 0., 0.) (0., 0., -1., 0., 0., 0.) (0., 0., 0., 0., -1., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0.) (1., 0., 0., 0., 0., 0.) (0., 0., 1., 0., 0., 0.) (0., 0., 0., 0., 1., 0.) (0., 1., 0., 0., 0., 0.) (0., 0., 0., 1., 0., 0.) (0., 0., 0., 0., 0., 1.);
        (1., 1., 0., 0., 0., 0.) (0., 0., 1., 1., 0., 0.) (0., 0., 0., 0., 1., 1.) (-3.0, -2., 0., 0., 0., 0.) (0., 0., -3.0, -2., 0., 0.) (0., 0., 0., 0., -3.0, -2.) (-1., -2., 0., 0., 0., 0.) (0., 0., -1., -2., 0., 0.) (0., 0., 0., 0., -1., -2.) (2., 0., 0., 0., 0., 0.) (0., 0., 2., 0., 0., 0.) (0., 0., 0., 0., 2., 0.) (1., 2., 0., 0., 0., 0.) (0., 0., 1., 2., 0., 0.) (0., 0., 0., 0., 1., 2.) (0., 1., 0., 0., 0., 0.) (0., 0., 0., 1., 0., 0.) (0., 0., 0., 0., 0., 1.);
        (-1., -1., 0., 0., 0., 0.) (0., 0., -1., -1., 0., 0.) (0., 0., 0., 0., -1., -1.) (0.6, -0.4, 0., 0., 0., 0.) (0., 0., 0.6, -0.4, 0., 0.) (0., 0., 0., 0., 0.6, -0.4) (-0.6, 0.4, 0., 0., 0., 0.) (0., 0., -0.6, 0.4, 0., 0.) (0., 0., 0., 0., -0.6, 0.4) (0.4, 0., 0., 0., 0., 0.) (0., 0., 0.4, 0., 0., 0.) (0., 0., 0., 0., 0.4, 0.) (0.6, 0.4, 0., 0., 0., 0.) (0., 0., 0.6, 0.4, 0., 0.) (0., 0., 0., 0., 0.6, 0.4) (0., 0.6, 0., 0., 0., 0.) (0., 0., 0., 0.6, 0., 0.) (0., 0., 0., 0., 0., 0.6)]
Hbx = H[(2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.) (-4.0, -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4.0, -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4.0, -2., -2., 0.) (0., -2., -2., -4.0, 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4.0, 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.0) (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.);
        (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.) (-4.0, -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4.0, -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4.0, -2., -2., 0.) (0., -2., -2., -4.0, 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4.0, 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.0) (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.);
        (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.) (-4.0, -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4.0, -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4.0, -2., -2., 0.) (0., -2., -2., -4.0, 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4.0, 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.0) (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.);
        (2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 2., 2., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2.) (-4.0, -2., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., -4.0, -2., -2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., -4.0, -2., -2., 0.) (0., -2., -2., -4.0, 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., -2., -2., -4.0, 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., -2., -2., -4.0) (2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0.) (0., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 2., 2., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 0.) (0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.) (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.)]
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])

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
Gbx = _Gbx(D,order,x,G)
Hbx = _Hbx(D,order,x,H)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=Gbx[1,:],gradgrad=Hbx[1,:])


############################################################################
# Tests for ND Bernstein polynomial with arbitrary barycentric coordinates #
############################################################################
#
T = Float64

D = 0
vertices = (Point(), )
b = BernsteinBasisOnSimplex(Val(D), Float64, 3, vertices)

D = 2

vertices = (Point(0.,1.), Point(0.,2.), Point(0.,3.))
@test_throws DomainError Polynomials._compute_cart_to_bary_matrix(vertices, Val(D+1))

vertices = (Point(5.,0.), Point(7.,2.), Point(0.,3.))

b = BernsteinBasisOnSimplex(Val(D), Float64, 3, vertices)
x = [Point(.0,.5), Point(1.,.5), Point(.2,.3), Point(5.,0.), Point(7.,2.), Point(0.,3.)]
x1 = x[1]

for xi in x
  λi = Polynomials._cart_to_bary(xi, b.cart_to_bary_matrix)
  @test sum(λi) ≈ 1.
  @test xi ≈ sum(λi .* vertices)
end


# Scalar value in 2D
D = 2
vertices = (Point(5.,0.), Point(7.,2.), Point(0.,3.))
x2λ = Polynomials._compute_cart_to_bary_matrix(vertices, Val(D+1))

G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 4
b = BernsteinBasisOnSimplex(Val(D),V,order,vertices)
bx  = _bx( D,order,x,  x2λ)
Gbx = _Gbx(D,order,x,G,x2λ)
Hbx = _Hbx(D,order,x,H,x2λ)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(b,x1,bx[1,:],≈,grad=Gbx[1,:],gradgrad=Hbx[1,:])


# scalar valued in 3D
D = 3
vertices = (Point(5.,0.,0.), Point(7.,0.,2.), Point(0.,3.,3.), Point(3.,0.,3.))
x2λ = Polynomials._compute_cart_to_bary_matrix(vertices, Val(D+1))

x = [Point(0.,0.,1.), Point(.5,.5,.5), Point(1.,.2,.4), Point(.2,.4,.3)]
x1 = x[1]

V = Float64
G = gradient_type(V,x1)
H = gradient_type(G,x1)

order = 4
b = BernsteinBasisOnSimplex(Val(D),V,order,vertices)
bx  = _bx( D,order,x,  x2λ)
Gbx = _Gbx(D,order,x,G,x2λ)
Hbx = _Hbx(D,order,x,H,x2λ)
test_field_array(b,x,bx,≈, grad=Gbx, gradgrad=Hbx)
test_field_array(b,x1,bx[1,:],≈,grad=Gbx[1,:],gradgrad=Hbx[1,:])

end # module
