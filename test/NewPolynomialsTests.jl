##
using Test
using Numa
using Numa.NewPolynomials
using Numa.FieldValues

using Numa.NewPolynomials: MultivariatePolynomialBasis
using Numa.NewPolynomials: UnivariatePolynomialBasis
using Numa.NewPolynomials: UnivariateMonomialBasis
using Numa.NewPolynomials: evaluate
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
using Numa.NewPolynomials: derivative
ders = derivative(a,points,numd=1)
@test ders[3,3][1] == 6.0
##
using Numa.NewPolynomials: gradient
grad = gradient(a)
length(grad)
# gradient(grad)
using Numa.NewPolynomials: evaluate
gval = evaluate(grad, points)
@test gval[3,3][1] == 6.0
##
using Numa.NewPolynomials: TensorProductMonomialBasis
D = 2
orders=[2,3]
tpmb = TensorProductMonomialBasis(orders)
length(tpmb)
@test typeof(tpmb) <: MultivariatePolynomialBasis{2,ScalarValue}
points = [ Point{2}(1.0, 1.0), Point{2}(1.0, 2.0), Point{2}(2.0, 2.0)]
vals = evaluate(tpmb, points)
@test vals[12,3] == 32.0
##
