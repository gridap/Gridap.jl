module IntegrateTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: OtherMockBasis, MockBasis, MockField
using FillArrays
using Gridap.TensorValues

p1 = Point(2,2)
p2 = Point(4,2)
p3 = Point(1,3)
p4 = Point(5,2)
x = [p1,p2,p3,p4]
np = length(x)

d = 2
v = 3.0
ndof = 8
i = MockField{d}(v)
w = fill(1/np,np)
j = TensorValue(1.0,0.0,0.0,1.0)

@test integrate(i,x,w,j) == 9.0

ri = 5.0
r = MockBasis{d}(ri,ndof)

@test integrate(r,x,w,j) == fill(15.0,ndof)

l = 10
ai = Fill(i,l)
aw = fill(w,l)
aj = fill(j,l)
ax = fill(x,l)
ar = Fill(r,l)

s = integrate(ai,ax,aw,aj)
test_array(s,fill(9.0,l))

s = integrate(ar,ax,aw,aj)
test_array(s,fill(fill(15.0,ndof),l))

c = fill(1.0,ndof)
f = OtherMockBasis{d}(ndof)
ϕ = lincomb(f,c)
j = ∇(ϕ)
w = fill(1/np,np)
b = attachmap(r,ϕ)

@test integrate(i,x,w,j) == 38016.0

@test integrate(b,x,w,j) == fill(63360.0,ndof)

l = 10
ac = fill(c,l)
af = Fill(f,l)
aϕ = lincomb(af,ac)
aj = ∇(aϕ)
aw = Fill(w,l)
ax = Fill(x,l)
ar = Fill(r,l)
ab = attachmap(ar,aϕ)

s = integrate(ab,ax,aw,aj)
test_array(s,fill(fill(63360.0,ndof),l))


end # module
