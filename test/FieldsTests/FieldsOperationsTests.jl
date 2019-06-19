module FieldsOperationsTests

using Test
using Gridap

using ..FieldsMocks

p1 = Point(1,2)
p2 = Point(0,3)
p3 = Point(7,4)
p4 = Point(8,1)
p = [p1,p2,p3,p4]

f = MockField(2,Int)
v = f + f

r = evaluate(f,p)
rv = r .+ r

r = evaluate(∇(f),p)
rg = r .+ r

test_field(v,p,rv,rg)

f = MockBasis(2,Int)
v = f + f

r = evaluate(f,p)
rv = r .+ r

r = evaluate(∇(f),p)
rg = r .+ r

test_basis(v,p,rv,rg)

change = ones(Int,3,3)
basis = MockBasis(2,Int)
v = change_basis(basis,change)

rv = change*evaluate(basis,p)
rg = change*evaluate(∇(basis),p)
test_basis(v,p,rv,rg)

end # module
