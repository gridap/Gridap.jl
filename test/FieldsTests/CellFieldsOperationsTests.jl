module CellFieldsOperationsTests

using Test
using Gridap
using ..CellFieldsMocks
using ..CellValuesMocks

import Gridap: ∇

l = 10
cf = IterCellFieldMock(2,Int,l)

p1 = Point(1,2)
p2 = Point(0,3)
p3 = Point(7,4)
p4 = Point(8,1)
p = [p1,p2,p3,p4]
cp = TestIterCellValue(p,l)

f,_ = iterate(cf)
fg = gradient(f)
r = evaluate(f,p)
rg = evaluate(fg,p)

cf2 = -cf
v = [-r for i in 1:l]
g = [-rg for i in 1:l]
test_iter_cell_field(cf2,cp,v,g)

cf = IndexCellFieldMock(2,Int,l)
cp = TestIndexCellValue(p,l)
cf2 = -cf
test_index_cell_field(cf2,cp,v,g)

cf = IterCellFieldMock(2,Int,l)
cf2 = cf+cf
v = [r+r for i in 1:l]
g = [rg+rg for i in 1:l]
test_iter_cell_field(cf2,cp,v,g)

cf = IndexCellFieldMock(2,Int,l)
cp = TestIndexCellValue(p,l)
cf2 = cf+cf
test_index_cell_field(cf2,cp,v,g)

cf = IterCellBasisMock(2,Int,l)
cp = TestIterCellValue(p,l)

f,_ = iterate(cf)
fg = gradient(f)
r = evaluate(f,p)
rg = evaluate(fg,p)
v = [r+r for i in 1:l]
g = [rg+rg for i in 1:l]
cf2 = cf+cf
test_iter_cell_basis(cf2,cp,v,g)

# varinner

cf = IterCellFieldMock(2,Int,l)
cb = IterCellBasisMock(2,Int,l)

f,_ = iterate(cf)
r = evaluate(f,p)
b,_ = iterate(cb)
rb = evaluate(b,p)

v = [ r.*r for i in 1:l ]
cs = varinner(cf,cf)
test_iter_cell_map(cs,cp,v)

ndofs, npoins = size(rb)
rv = zeros(eltype(r),ndofs,npoins)
for j in 1:npoins
  for i in 1:ndofs
    rv[i,j] = rb[i,j]*r[j]
  end
end
v = [ rv for i in 1:l ]

cs = varinner(cb,cf)
test_iter_cell_map(cs,cp,v)

rv = zeros(eltype(r),ndofs,ndofs,npoins)
for k in 1:npoins
  for j in 1:ndofs
    for i in 1:ndofs
      rv[i,j,k] = rb[i,k]*rb[j,k]
    end
  end
end
v = [ rv for i in 1:l ]

cs = varinner(cb,cb)
test_iter_cell_map(cs,cp,v)

# lincomb

cf = IterCellFieldMock(2,Float64,l)
cb = IterCellBasisMock(2,Float64,l)

f,_ = iterate(cf)
r = evaluate(f,p)
b,_ = iterate(cb)
rb = evaluate(b,p)

u = rand(eltype(r),ndofs)
cu = TestIterCellValue(u,l)

ndofs, npoins = size(rb)
rv = zeros(eltype(r),npoins)
for j in 1:npoins
  for i in 1:ndofs
    rv[j] += rb[i,j]*u[i]
  end
end
v = [ rv for i in 1:l ]

rbg = evaluate(∇(b),p)
rv = zeros(eltype(rbg),npoins)
for j in 1:npoins
  for i in 1:ndofs
    rv[j] += rbg[i,j]*u[i]
  end
end
g = [ rv for i in 1:l ]

cs = lincomb(cb,cu)
test_iter_cell_field(cs,cp,v,g)

# compose

cf = IterCellFieldMock(2,Int,l)
ufun(x) = 2*x[1]

f,_ = iterate(cf)
r = evaluate(f,p)
v = [ ufun.(r) for i in 1:l ]

cr = compose(ufun,cf)
test_iter_cell_field_without_grad(cr,cp,v)

ufungrad(x) = VectorValue(2,0)
∇(::typeof(ufun)) = ufungrad

vg = [ ufungrad.(r) for i in 1:l ]

cr = compose(ufun,cf)
test_iter_cell_field(cr,cp,v,vg)

# attachgeomap

cf = IterCellGeomapMock(2,Float64,l)
cb = IterCellBasisMock(2,Float64,l)

b,_ = iterate(cb)
rb = evaluate(b,p)
rbg = evaluate(∇(b),p)

f,_ = iterate(cf)
r = evaluate(∇(f),p)

ndofs, npoins = size(rb)
v = [ rb for i in 1:l ]
g = [ reshape(inv.(r),1,npoins) .* rbg for i in 1:l ]

cb2 = attachgeomap(cb,cf)

test_iter_cell_basis(cb2,cp,v,g)

end # module

