module CellArrayApplyTests

using Test
using Gridap
using Gridap.CachedArrays
using ..CellValuesMocks
using TensorValues

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

end # module
