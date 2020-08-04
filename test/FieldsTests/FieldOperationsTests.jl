module FieldOperationsTests

using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: MockField, MockBasis
using Gridap.TensorValues
using FillArrays

p1 = Point(2,2)
p2 = Point(4,2)
p3 = Point(1,3)
p4 = Point(5,2)
x = [p1,p2,p3,p4]
np = length(x)

va = 3.0
vb = 3.3
d = 2
fa = MockField{d}(va)
fb = MockField{d}(vb)
fax = evaluate(fa,x)
fbx = evaluate(fb,x)
∇fax = evaluate(gradient(fa),x)
∇fbx = evaluate(gradient(fb),x)

wa = 4.5
wb = 4.5
ndofa = 8
ndofb = 5
ba = MockBasis{d}(wa,ndofa)
bb = MockBasis{d}(wb,ndofb)
bax = evaluate(ba,x)
bbx = evaluate(bb,x)
∇bax = evaluate(gradient(ba),x)
∇bbx = evaluate(gradient(bb),x)

l = 10
ax = fill(x,l)
afa = Fill(fa,l)
afb = Fill(fb,l)
aba = Fill(ba,l)
abb = Fill(bb,l)

g = field_operation(inner,fa,fb)
gx = [inner(ai,bj) for (ai,bj) in zip(fax,fbx)]
test_field(g,x,gx)

ag = field_array_operation(inner,afa,afb)
agx = fill(gx,l)
test_array_of_fields(ag,ax,agx)

struct FieldPlaceHolder <: Field end

ag = field_array_operation(FieldPlaceHolder,inner,fill(fa,l),afb)
test_array(evaluate_field_array(ag,ax),agx)

g = field_operation(inner,ba,fb)
gx = zeros(size(bax))
for p in 1:np
  for i in 1:ndofa
    gx[p,i] = inner(bax[p,i],fbx[p])
  end
end
test_field(g,x,gx)

ag = field_array_operation(inner,aba,afb)
agx = fill(gx,l)
test_array_of_fields(ag,ax,agx)

g = field_operation(inner,fa,bb)
gx = zeros(size(bbx))
for p in 1:np
  for i in 1:ndofb
    gx[p,i] = inner(fax[p],bbx[p,i])
  end
end
test_field(g,x,gx)

ag = field_array_operation(inner,afa,abb)
agx = fill(gx,l)
test_array_of_fields(ag,ax,agx)

g = field_operation(inner,ba,bb)
gx = zeros(np,ndofa,ndofb)
for p in 1:np
  for i in 1:ndofa
    for j in 1:ndofb
    gx[p,i,j] = inner(bax[p,i],bbx[p,j])
    end
  end
end
test_field(g,x,gx)

ag = field_array_operation(inner,aba,abb)
agx = fill(gx,l)
test_array_of_fields(ag,ax,agx)

g = field_operation(+,fa,fb)
gx = [+(ai,bj) for (ai,bj) in zip(fax,fbx)]
∇gx = [+(ai,bj) for (ai,bj) in zip(∇fax,∇fbx)]
test_field(g,x,gx,grad=∇gx)

ag = field_array_operation(+,afa,afb)
agx = fill(gx,l)
∇agx = fill(∇gx,l)
test_array_of_fields(ag,ax,agx,grad=∇agx)

g = field_operation(+,fa)
gx = [+ai for ai in fax]
∇gx = [+ai for ai in ∇fax]
test_field(g,x,gx,grad=∇gx)

ag = field_array_operation(+,afa)
agx = fill(gx,l)
∇agx = fill(∇gx,l)
test_array_of_fields(ag,ax,agx,grad=∇agx)

g = field_operation(*,va,fa)
gx = [*(va,ai) for ai in fax]
∇gx = [*(va,ai) for ai in ∇fax]
test_field(g,x,gx,grad=∇gx)


end # module
