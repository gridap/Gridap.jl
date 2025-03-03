module ChangeBasisTests

using Test
using Gridap.Fields
using Gridap.Polynomials
using StaticArrays

xi = Point(0.)
np = 5
x = fill(xi,np)

######################
# Polynomial interface
######################

struct MockPolynomial <: Polynomial end
@test_throws ErrorException isHierarchical(Polynomial)

# Interfaces to implement
@test_throws ErrorException isHierarchical(MockPolynomial)

D = 1
K = 0
c = zero(MMatrix{D,K+1})

@test_throws ErrorException Polynomials._evaluate_1d!(MockPolynomial, Val(1), c, xi, 1)
@test_throws ErrorException Polynomials._gradient_1d!(MockPolynomial, Val(1), c, xi, 1)
@test_throws ErrorException Polynomials._hessian_1d!( MockPolynomial, Val(1), c, xi, 1)
@test_throws ErrorException Polynomials._derivatives_1d!(MockPolynomial, Val(1), (nothing,nothing,nothing,nothing), xi, 1)

function Polynomials._evaluate_1d!(
  ::Type{MockPolynomial},::Val{K}, cc::AbstractMatrix{T}, xi, d) where {K,T<:Number}

  cc[1,1] = 1.
end
Polynomials._derivatives_1d!(MockPolynomial, Val(1), (c,), xi, 1)
@test c[1][1] == 1.

###########################
# PolynomialBasis interface
###########################

T = Float64
D = 1
struct MockPolyBasis <: PolynomialBasis{D,T,0,MockPolynomial} end

mb = MockPolyBasis()


# Implemented interfaces
@test IndexStyle(mb) == IndexLinear()
@test return_type(mb) == T
@test get_order(mb) == 0
@test mb[1] == MockPolynomial()


# Interfaces to implement
@test_throws ErrorException size(mb)
import Base.size
Base.size(::MockPolyBasis) = (1,)
@test length(mb) == 1


r, c = return_cache(mb,x)

@test_throws ErrorException Polynomials._evaluate_nd!(mb, xi, r, 1, c)

∇mb = FieldGradientArray{1}(mb)
r, c, g = return_cache(∇mb,x)
s = MVector{D,T}(0.)

@test_throws ErrorException Polynomials._gradient_nd!(mb, xi, r, 1, c, g, s)

Hmb = FieldGradientArray{2}(mb)
r, c, g, h = return_cache(Hmb,x)
s = MMatrix{D,D,T}(0.)

@test_throws ErrorException Polynomials._hessian_nd!(mb, xi, r, 1, c, g, h, s)


end
