module ExteriorCalculusBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Polynomials: _p_filter, _q_filter, _ser_filter
using ForwardDiff
using StaticArrays

# Source: https://www-users.cse.umn.edu/~arnold/femtable/background.html
FEEC_length(r,k,D,::Val{:P⁻}) = binomial(r+D,r+k)*binomial(r+k-1,k)
FEEC_length(r,k,D,::Val{:P }) = binomial(r+D,r+k)*binomial(r+k,k)
FEEC_length(r,k,D,::Val{:Q⁻}) = binomial(D,k) * r^k * (r+1)^(D-k)
FEEC_length(r,k,D,::Val{:S }) = sum( 2^(D-d)*binomial(D,d)*binomial(r-d+2k,d)*binomial(d,k) for d in k:D )

T = Float64   # scalar type
r = 3         # Polynomial order

no_hessian = true
@inline function _test_bases(b::FEECPolyBasis, b2, no_hessian=false)
  #D = get_dimension(b)
  r = get_FEEC_poly_degree(b)
  k = get_FEEC_form_degree(b)
  F = get_FEEC_family(b)
  @test length(b) == FEEC_length(r,k,D,Val(F))
  @test b._basis isa typeof(b2)
  @test evaluate(b,x) == evaluate(b2,x)
  @test evaluate(Broadcasting(∇ )(b),x) == evaluate(Broadcasting(∇ )(b2),x)
  if no_hessian
    @test_throws ErrorException evaluate(Broadcasting(∇∇)(b),x) == evaluate(Broadcasting(∇∇)(b2),x)
  else
    @test evaluate(Broadcasting(∇∇)(b),x) == evaluate(Broadcasting(∇∇)(b2),x)
  end
end

# 1D
D = 1
V = T
x = [Point(0.),Point(1.),Point(.4)]

k = 0
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2)


k = 1
V = T
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_ser_filter)
_test_bases(b,b2)


# 2D
D = 2
V = T
x = [Point(0.,0.),Point(1.,0),Point(0,.4),Point(.6,.4)]

k = 0
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2)


k = 1
V = VectorValue{D,T}
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = PCurlGradBasis(Monomial,Val(D),T,r-1)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,no_hessian)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = QCurlGradBasis(Monomial,Val(D),T,r-1)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2)

@test_throws ErrorException FEECPolyBasis(Val(D),T,r,k,:S)
# b = FEECPolyBasis(Val(D),T,r,k,:S)
# b2 = # TODO
# @test b._basis isa PolynomialBasis{D,V,r,Monomial}
#_test_bases(b,b2)


k = 2
V = T
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2)


# 3D
D = 3
V = T
x = [Point(0.,0,0),Point(1.,0,0),Point(0,.4,.3),Point(.6,.4,.5)]

k = 0
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2)


k = 1
V = VectorValue{D,T}
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = PGradBasis(Monomial,Val(D),T,r-1)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,no_hessian)

@test_throws ErrorException FEECPolyBasis(Val(D),T,r,k,:P)
# b = FEECPolyBasis(Val(D),T,r,k,:P)
# b2 =
# @test b._basis isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = QGradBasis(Monomial,Val(D),T,r-1)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2)

@test_throws ErrorException FEECPolyBasis(Val(D),T,r,k,:S)
# b = FEECPolyBasis(Val(D),T,r,k,:P)
# b2 =
# @test b._basis isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2)


k = 2
V = VectorValue{D,T}
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = PCurlGradBasis(Monomial,Val(D),T,r-1)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,no_hessian)

@test_throws ErrorException FEECPolyBasis(Val(D),T,r,k,:P)
# b = FEECPolyBasis(Val(D),T,r,k,:P)
# b2 =
# @test b._basis isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = QCurlGradBasis(Monomial,Val(D),T,r-1)
@test b._basis isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2)

@test_throws ErrorException FEECPolyBasis(Val(D),T,r,k,:S)
# b = FEECPolyBasis(Val(D),T,r,k,:P)
# b2 =
# @test b._basis isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2)


k = 3
V = T
b = FEECPolyBasis(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2)

b = FEECPolyBasis(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2)


end # module
