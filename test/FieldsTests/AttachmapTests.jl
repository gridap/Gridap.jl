module AttachmapTests

using Test
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: OtherMockBasis, MockBasis
using FillArrays
using LinearAlgebra

p1 = Point(2,2)
p2 = Point(4,2)
p3 = Point(1,3)
p4 = Point(5,2)
x = [p1,p2,p3,p4]

np = length(x)
v = 2.0
d = 2
ndof = 8
wi = 3.0
w = fill(wi,ndof)
f = OtherMockBasis{d}(ndof)
r = MockBasis{d}(v,ndof)
ϕ = lincomb(f,w)

b = attachmap(r,ϕ)
rx = evaluate(r,x)
bx = rx
∇rx = evaluate(∇(r),x)
∇ϕx = evaluate(∇(ϕ),x)
∇bx = zeros(VectorValue{d,Float64},np,ndof)
for i in 1:np
  jacinv = inv(∇ϕx[i])
  for j in 1:ndof
    ∇bx[i,j] = jacinv⋅∇rx[i,j]
  end
end
test_field(b,x,bx,grad=∇bx)

l = 10
ax = fill(x,l)
af = Fill(f,l)
ar = Fill(r,l)
aw = fill(w,l)
aϕ = lincomb(af,aw)

ab = attachmap(ar,aϕ)
abx = fill(bx,l)
a∇bx = fill(∇bx,l)
test_array_of_fields(ab,ax,abx,grad=a∇bx)

end # module
