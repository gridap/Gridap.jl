##
using Test
using Numa
using Numa.Polynomials
using Numa.FieldValues

using Numa.Polynomials: MultivariatePolynomialBasis
using Numa.Polynomials: UnivariatePolynomialBasis
using Numa.Polynomials: UnivariateMonomialBasis
using Numa.Polynomials: evaluate
##

@test UnivariatePolynomialBasis <: MultivariatePolynomialBasis{1,ScalarValue}
a = UnivariateMonomialBasis(2)
@test typeof(a) <: MultivariatePolynomialBasis{1,ScalarValue}
@test length(a) == 3
#evaluate UnivariateMonomialBasis
points = [Point{1}(1.0), Point{1}(2.0), Point{1}(3.0)]
values = evaluate(a,points)
@test values[3,3] == 9
##
using Numa.Polynomials: derivative
ders = derivative(a,points,numd=1)
@test ders[3,3][1] == 6.0
##
using Numa.Polynomials: gradient
grad = gradient(a)
length(grad)
# gradient(grad)
using Numa.Polynomials: evaluate
gval = evaluate(grad, points)
@test gval[3,3][1] == 6.0
##
using Numa.Polynomials: TensorProductMonomialBasis
D = 2
orders=[2,3]
tpmb = TensorProductMonomialBasis{D,ScalarValue}(orders)
@test typeof(tpmb) <: MultivariatePolynomialBasis{2,ScalarValue}
points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
vals = evaluate(tpmb, points)
@test vals[12,3] == 32.0
##

##
using Numa.Polynomials: TensorProductMonomialBasis
D = 2
orders=[2,3]
tpmb = TensorProductMonomialBasis{D,VectorValue{D}}(orders)
@test typeof(tpmb) <: MultivariatePolynomialBasis{2,VectorValue{D}}
points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
vals = evaluate(tpmb, points)
@test vals[12,3][1] == 32.0
@test vals[12,3][2] == 0.0
@test vals[24,3][1] == 0.0
@test vals[24,3][2] == 32.0
##

##
using Numa.Polynomials: TensorProductMonomialBasis
D = 2
orders=[2,3]
tpmb = TensorProductMonomialBasis{D,TensorValue{D}}(orders)
@test typeof(tpmb) <: MultivariatePolynomialBasis{2,TensorValue{D}}
points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
vals = evaluate(tpmb, points)
@test vals[12,3][1,1] == 32.0
@test vals[12,3][2,1] == 0.0
@test vals[12,3][1,2] == 0.0
@test vals[12,3][2,2] == 0.0
@test vals[24,3][1,1] == 0.0
@test vals[24,3][2,1] == 32.0
@test vals[24,3][1,2] == 0.0
@test vals[24,3][2,2] == 0.0
@test vals[36,3][1,1] == 0.0
@test vals[36,3][2,1] == 0.0
@test vals[36,3][1,2] == 32.0
@test vals[36,3][2,2] == 0.0
@test vals[48,3][1,1] == 0.0
@test vals[48,3][2,1] == 0.0
@test vals[48,3][1,2] == 0.0
@test vals[48,3][2,2] == 32.0
##
