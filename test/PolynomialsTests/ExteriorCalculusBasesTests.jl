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

@inline function _test_bases(b, b2, r,k,F,D, no_hessian=false)
  @test length(b) == FEEC_length(r,k,D,Val(F))
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
V = T
x = [Point(0.),Point(1.),Point(.4)]

D, k = 1, 0
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2,r,k,:S,D)


D, k = 1, 1
V = T
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_ser_filter)
_test_bases(b,b2,r,k,:S,D)


# 2D
V = T
x = [Point(0.,0.),Point(1.,0),Point(0,.4),Point(.6,.4)]

D, k = 2, 0
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P ,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2,r,k,:S ,D)


D, k = 2, 1
V = VectorValue{D,T}
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = PCurlGradBasis(Monomial,Val(D),T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = QCurlGradBasis(Monomial,Val(D),T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
# b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
# b2 = # TODO
# @test b isa PolynomialBasis{D,V,r,Monomial}
#_test_bases(b,b3,r,k,:S,D)


D, k = 2, 2
V = T
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:S,D)


# 3D
V = T
x = [Point(0.,0,0),Point(1.,0,0),Point(0,.4,.3),Point(.6,.4,.5)]

D, k = 3, 0
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2,r,k,:S,D)


D, k = 3, 1
V = VectorValue{D,T}
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = PGradBasis(Monomial,Val(D),T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

@test_throws ErrorException FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
# b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
# b2 =
# @test b isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = QGradBasis(Monomial,Val(D),T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
# b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
# b2 =
# @test b isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2,r,k,:S,D)


D, k = 3, 2
V = VectorValue{D,T}
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = PCurlGradBasis(Monomial,Val(D),T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:P⁻,D,no_hessian)

@test_throws ErrorException FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
# b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
# b2 =
# @test b isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = QCurlGradBasis(Monomial,Val(D),T,r-1)
@test b isa PolynomialBasis{D,V,Monomial}
_test_bases(b,b2,r,k,:Q⁻,D)

@test_throws ErrorException FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
# b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
# b2 =
# @test b isa PolynomialBasis{D,V,r,Monomial}
# _test_bases(b,b2,r,k,:S,D)


D, k = 3, 3
V = T
b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:P)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:Q⁻)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEECPolyBasis_trampoline(Val(D),T,r,k,:S)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:S,D)


end # module
