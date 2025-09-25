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
@test Polynomials._default_poly_type(:Q⁻) == Legendre
@test Polynomials._default_poly_type(:S)  == Legendre
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
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2,r,k,:S,D)


D, k = 1, 1
V = T
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_ser_filter)
_test_bases(b,b2,r,k,:S,D)


# 2D
V = T
x = [Point(0.,0.),Point(1.,0),Point(0,.4),Point(.6,.4)]

D, k = 2, 0
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P ,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2,r,k,:S ,D)


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
V = T
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:S,D)


# 3D
V = T
x = [Point(0.,0,0),Point(1.,0,0),Point(0,.4,.3),Point(.6,.4,.5)]

D, k = 3, 0
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),V,r,_ser_filter)
_test_bases(b,b2,r,k,:S,D)


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
V = T
b = FEEC_poly_basis(Val(D),T,r,k,:P⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_p_filter)
_test_bases(b,b2,r,k,:P⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:P,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:P,D)

b = FEEC_poly_basis(Val(D),T,r,k,:Q⁻,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r-1,_q_filter)
_test_bases(b,b2,r,k,:Q⁻,D)

b = FEEC_poly_basis(Val(D),T,r,k,:S,Monomial)
b2 = CartProdPolyBasis(Monomial,Val(D),T,r,_p_filter)
_test_bases(b,b2,r,k,:S,D)


@test_throws ErrorException FEEC_poly_basis(Val(4),T,r,1,:P⁻,Monomial)
@test_throws ErrorException FEEC_poly_basis(Val(4),T,r,2,:P⁻,Monomial)

end # module
