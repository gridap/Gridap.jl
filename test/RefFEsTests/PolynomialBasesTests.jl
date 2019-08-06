module PolynomialBasesTests

using Test
using Gridap

p1 = Point(1,2)
p2 = Point(0,3)
p = [p1,p2]

T = VectorValue{3,Float64}
orders = (1,1)
basis = MonomialBasis(T,orders)

v = evaluate(basis,p)
g = evaluate(∇(basis),p)

test_basis(basis,p,v,g)

T = Float64
orders = (1,1)
basis = MonomialBasis(T,orders)

v = evaluate(basis,p)
g = evaluate(∇(basis),p)

test_basis(basis,p,v,g)

T = VectorValue{2,Float64}
filter(e,order) = sum(e) <= order
order = 1
basis = MonomialBasis(2,T,filter,order)

@test length(basis) == 6

T = VectorValue{2,Float64}
basis = GradMonomialBasis(T,2)
v = evaluate(basis,p)
g = evaluate(∇(basis),p)

@test length(basis) == 12

p = Point{3,Int}[(1,2,3),(3,2,4)]

T = VectorValue{3,Float64}
order = 1
basis = CurlGradMonomialBasis(T,order)
v = evaluate(basis,p)
g = evaluate(∇(basis),p)

@test length(basis) == 6

end # module
