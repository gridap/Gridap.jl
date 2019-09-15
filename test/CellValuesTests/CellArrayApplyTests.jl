module CellArrayApplyTests

using Test
using Gridap
using Gridap.CachedArrays
using ..CellValuesMocks
using TensorValues
using Gridap.Kernels: IntegrateNumberKernel
using Gridap.Kernels: IntegrateArrayKernel

l = 10

a = [1,2,3]
v = TestIterCellValue(a,l)
w = apply(-,v,broadcast=true)
o = [ CachedArray(-vi) for vi in v ] 
test_iter_cell_array(w,o)

oo = collect(w)
@test oo == o

a = [1,2,3]
b = [3,2,1]
u = TestIterCellValue(a,l)
v = TestIterCellValue(b,l)
w = apply(-,u,v,broadcast=true)
o = [ CachedArray(ui.-vi) for (ui,vi) in zip(u,v) ] 
test_iter_cell_array(w,o)

v1 = VectorValue(2,3)
v2 = VectorValue(3,2)
v3 = VectorValue(1,2)

ax = [rand(1,3,4), 1.0      , [v1,v2,v3], [v1,v2,v3], 1         ]
bx = [rand(2,3,1), rand(2,3), [v2,v3,v1], v1        , [v2,v3,v1]]

for (a,b) in zip(ax,bx)
  u = TestIterCellValue(a,l)
  v = TestIterCellValue(b,l)
  o = [ CachedArray(ui.-vi) for (ui,vi) in zip(u,v) ] 
  w = apply(-,u,v,broadcast=true)
  test_iter_cell_array(w,o)
end

a = [1,2,3]
v = TestIndexCellValue(a,l)
w = apply(-,v,broadcast=true)
o = [ CachedArray(-vi) for vi in v ] 
test_index_cell_array(w,o)

a = [1,2,3]
b = [3,2,1]
u = TestIndexCellValue(a,l)
v = TestIndexCellValue(b,l)
w = apply(-,u,v,broadcast=true)
o = [ CachedArray(ui.-vi) for (ui,vi) in zip(u,v) ] 
test_index_cell_array(w,o)

v1 = VectorValue(2,3)
v2 = VectorValue(3,2)
v3 = VectorValue(1,2)

ax = [rand(1,3,4), 1.0      , [v1,v2,v3], [v1,v2,v3], 1         ]
bx = [rand(2,3,1), rand(2,3), [v2,v3,v1], v1        , [v2,v3,v1]]

for (a,b) in zip(ax,bx)
  u = TestIndexCellValue(a,l)
  v = TestIndexCellValue(b,l)
  o = [ CachedArray(ui.-vi) for (ui,vi) in zip(u,v) ] 
  w = apply(-,u,v,broadcast=true)
  test_index_cell_array(w,o)
end

a = [1,2,3]
b = [3,2,1]
c = [3,2,4]
u = TestIndexCellValue(a,l)
v = TestIndexCellValue(b,l)
w = TestIndexCellValue(c,l)
k = IntegrateNumberKernel{Int}()
z = apply(k,u,v,w)
r = fill(sum(a.*b.*c),l)
test_index_cell_number(z,r)

a = rand(3,2,4)
b = rand(4)
c = rand(4)
u = TestIndexCellValue(a,l)
v = TestIndexCellValue(b,l)
w = TestIndexCellValue(c,l)
k = IntegrateArrayKernel()
z = apply(k,u,v,w)
ri =  reshape( sum(a .* reshape(b,(1,1,4)) .* reshape(c,(1,1,4)), dims=3), (3,2))
r = fill(collect(ri),l)
test_index_cell_array(z,r)

end # module
