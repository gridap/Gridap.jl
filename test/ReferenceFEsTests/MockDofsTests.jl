module MockDofsTests

using Test
using Gridap.Fields
using Gridap.Fields: MockField
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: MockDofBasis

np = 4
p = Point(1,2)
x = fill(p,np)

v = 3.0
d = 2
f = MockField{d}(v)
fx = evaluate(f,x)

b = MockDofBasis(x)
test_dof(b,f,fx)

l = 10
ab = fill(b,l)
af = fill(f,l)
ax = fill(x,l)

afx = evaluate(af,ax)
abf = evaluate(ab,af)
@test abf == afx

end # module
