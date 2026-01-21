module ExteriorCalculusBasesTests

using Base: diff_names
using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Polynomials: _p_filter, _q_filter, _ser_filter
using ForwardDiff
using StaticArrays

@test Polynomials._default_poly_type(:P⁻) == Bernstein
@test Polynomials._default_poly_type(:P)  == Bernstein
@test Polynomials._default_poly_type(:Q⁻) == Polynomials.ModalC0
@test Polynomials._default_poly_type(:S)  == Polynomials.ModalC0
@test Polynomials._default_poly_type(:default) == Monomial

# differential geometry / exterior calculus isn't implemented yet
DG_calc = true
@test_throws ErrorException FEEC_space_definition_checks(Val(0),Float64,0,0,:P,false,DG_calc)

# no vector proxy for 1 < k < D-1 forms
D,k = 4, 2
@test_throws ErrorException FEEC_space_definition_checks(Val(D),Float64,0,k,:P)

# rotate_90 is for 2D 1-forms
rotate_90 = true
@test_warn "`rotate_90` kwarg" FEEC_space_definition_checks(Val(2),Float64,0,2,:P,rotate_90)
@test_warn "`rotate_90` kwarg" FEEC_space_definition_checks(Val(3),Float64,0,1,:P,rotate_90)

# cart_prod only for 0/D forms
cart_prod = true
rotate_90 = false
@test_throws "Cartesian product" FEEC_space_definition_checks(Val(2),Float64,0,1,:P,rotate_90; cart_prod)
@test FEEC_space_definition_checks(Val(2),VectorValue{3,Float64},0,0,:P,rotate_90; cart_prod)

# Source: https://www-users.cse.umn.edu/~arnold/femtable/background.html
FEEC_length(r,k,D,::Val{:P⁻}) = binomial(r+D,r+k)*binomial(r+k-1,k)
FEEC_length(r,k,D,::Val{:P }) = binomial(r+D,r+k)*binomial(r+k,k)
FEEC_length(r,k,D,::Val{:Q⁻}) = binomial(D,k) * r^k * (r+1)^(D-k)
FEEC_length(r,k,D,::Val{:S }) = sum( 2^(D-d)*binomial(D,d)*binomial(r-d+2k,d)*binomial(d,k) for d in k:D )

T = Float64   # scalar type
r = 3         # Polynomial order
no_hessian = true

@noinline function _test_bases(b, b2, r,k,F,D, no_hessian=false; cart_prod=false)
  ncomp = cart_prod ? num_indep_components(value_type(b)) : 1
  @test length(b) == FEEC_length(r,k,D,Val(F))*ncomp
  @test b isa typeof(b2)
  @test evaluate(b,x) == evaluate(b2,x)
  @test evaluate(Broadcasting(∇)(b),x) == evaluate(Broadcasting(∇)(b2),x)
  @test testvalue(typeof(b)) isa typeof(b)
  if no_hessian
    @test_throws ErrorException evaluate(Broadcasting(∇∇)(b),x) == evaluate(Broadcasting(∇∇)(b2),x)
  else
    @test evaluate(Broadcasting(∇∇)(b),x) == evaluate(Broadcasting(∇∇)(b2),x)
  end
end


# 1D
x = [Point(0.),Point(1.),Point(.4)]

D, k = 1, 0
function _k0D_tests(D,k,r,T)
  VC = VectorValue{2,T}
  m1     = k==D ? 1 : 0
  filter = k==D ? _p_filter : _ser_filter
  for Vi in (T, VC)
    cart_prod = Vi <: MultiValue

    # Monomial
    b = FEEC_poly_basis(Val(D),Vi,r,k,:P⁻,Monomial; cart_prod)
    b2 = CartProdPolyBasis(Monomial,Val(D),Vi,r-m1,_p_filter)
    _test_bases(b,b2,r,k,:P⁻,D; cart_prod)

    b = FEEC_poly_basis(Val(D),Vi,r,k,:P,Monomial; cart_prod)
    b2 = CartProdPolyBasis(Monomial,Val(D),Vi,r,_p_filter)
    _test_bases(b,b2,r,k,:P,D; cart_prod)

    b = FEEC_poly_basis(Val(D),Vi,r,k,:Q⁻,Monomial; cart_prod)
    b2 = CartProdPolyBasis(Monomial,Val(D),Vi,r-m1,_q_filter)
    _test_bases(b,b2,r,k,:Q⁻,D; cart_prod)

    b = FEEC_poly_basis(Val(D),Vi,r,k,:S,Monomial; cart_prod)
    b2 = CartProdPolyBasis(Monomial,Val(D),Vi,r,filter)
    _test_bases(b,b2,r,k,:S,D; cart_prod)

    # Bernstein not PmΛ
    if cart_prod
    b = FEEC_poly_basis(Val(D),Vi,r,k,:P⁻,Bernstein; cart_prod)
    b2 = BernsteinBasisOnSimplex{D}(Vi,r-m1)
    _test_bases(b,b2,r,k,:P⁻,D; cart_prod)

    b = FEEC_poly_basis(Val(D),Vi,r,k,:P,Bernstein; cart_prod)
    b2 = BernsteinBasisOnSimplex{D}(Vi,r)
    _test_bases(b,b2,r,k,:P,D; cart_prod)

    b = FEEC_poly_basis(Val(D),Vi,r,k,:Q⁻,Bernstein; cart_prod)
    b2 = CartProdPolyBasis(Bernstein,Val(D),Vi,r-m1,_q_filter)
    _test_bases(b,b2,r,k,:Q⁻,D; cart_prod)
    end

    if k==0
      @test_throws "hierarchical" FEEC_poly_basis(Val(D),Vi,r,k,:S,Bernstein; cart_prod)
    else
      b = FEEC_poly_basis(Val(D),Vi,r,k,:S,Bernstein; cart_prod)
      b2 = BernsteinBasisOnSimplex{D}(Vi,r)
      _test_bases(b,b2,r,k,:S,D; cart_prod)
    end
  end
end
_k0D_tests(D,k,r,T)

D, k = 1, 1
_k0D_tests(D,k,r,T)


# 2D
x = [Point(0.,0.),Point(1.,0),Point(0,.4),Point(.6,.4)]

D, k = 2, 0
_k0D_tests(D,k,r,T)


D, k = 2, 1
V = VectorValue{D,T}

b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = NedelecPolyBasisOnSimplex{D}(Monomial,T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
orders = [ r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D ]
b2 = CompWiseTensorPolyBasis{D}(Monomial, V, orders)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
# b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
# b2 = # TODO
# @test b isa PolynomialBasis{D,V,r,Monomial}
#_test_bases(b,b3,r,k,:S,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial; rotate_90=true)
b2 = RaviartThomasPolyBasis{D}(Monomial,T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial; rotate_90=true)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial; rotate_90=true)
orders = [ r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D ]
b2 = CompWiseTensorPolyBasis{D}(Monomial,V,orders)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEEC_poly_basis(Val(D),T,r,k,:S,Monomial; rotate_90=true)
# b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
# b2 = # TODO
# @test b isa PolynomialBasis{D,V,r,Monomial}
#_test_bases(b,b3,r,k,:S,D)


D, k = 2, 2
_k0D_tests(D,k,r,T)


# 3D
V = T
x = [Point(0.,0,0),Point(1.,0,0),Point(0,.4,.3),Point(.6,.4,.5)]

D, k = 3, 0
_k0D_tests(D,k,r,T)

D, k = 3, 1
V = VectorValue{D,T}
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = NedelecPolyBasisOnSimplex{D}(Monomial,T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = MonomialBasis(Val(D),V,r,Polynomials._p_filter)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
orders = [ r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D ]
b2 = CompWiseTensorPolyBasis{D}(Monomial, V, orders)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
# b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
# b2 =
# @test b isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2,r,k,:S,D)


D, k = 3, 2
V = VectorValue{D,T}
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = RaviartThomasPolyBasis{D}(Monomial,T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = MonomialBasis(Val(D),V,r,Polynomials._p_filter)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
orders = [ r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D ]
b2 = CompWiseTensorPolyBasis{D}(Monomial, V, orders)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
# b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
# b2 =
# @test b isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2,r,k,:S,D)


D, k = 3, 3
_k0D_tests(D,k,r,T)

@test_throws ErrorException FEEC_poly_basis(Val(4),T,r,1,:P⁻,Monomial)
@test_throws ErrorException FEEC_poly_basis(Val(4),T,r,2,:P⁻,Monomial)

end # module
