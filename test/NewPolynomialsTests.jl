##
using Test
using Numa
using Numa.NewPolynomials
using Numa.FieldValues

using Numa.NewPolynomials: MultivariatePolynomialBasis
using Numa.NewPolynomials: UnivariatePolynomialBasis
using Numa.NewPolynomials: UnivariateMonomialBasis
##

@test UnivariatePolynomialBasis <: MultivariatePolynomialBasis{1,ScalarValue}
a = UnivariateMonomialBasis(2)
@test typeof(a) <: MultivariatePolynomialBasis{1,ScalarValue}
@test length(a) == 3
#evaluate UnivariateMonomialBasis
points = [Point{1}(1.0), Point{1}(2.0), Point{1}(3.0)]
values = a(points)
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














# l=10000000
# println("+++ CellFieldsBench ( length = $l ) +++")
# print("Array of 1D points ->"); @time dopointmon(l)
# print("Array of Ints ->"); @time doarraymon(l)
