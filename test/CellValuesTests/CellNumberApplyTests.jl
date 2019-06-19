module CellNumberApplyTests

using Test
using Gridap
using ..CellValuesMocks
using StaticArrays
using TensorValues

l = 10

ax = [1.2, SVector(1.3,2.0), VectorValue(1.3,2.0), VectorValue(1.3,2.0)]
bx = [1, SVector(1.3,2.3), VectorValue(1.1,2.2),1]

for a in ax
  v = TestIterCellValue(a,l)
  w = apply(-,v)
  o = [ -vi for vi in v ] 
  test_iter_cell_value(w,o)
end

for (a,b) in zip(ax,bx)
  u = TestIterCellValue(a,l)
  v = TestIterCellValue(b,l)
  w = apply(-,u,v)
  o = [ ui-vi for (ui,vi) in zip(u,v) ] 
  test_iter_cell_value(w,o)
end

for a in [[1,2,3], SVector(1.3,2.0), VectorValue(1.3,2.0)]
  v = TestIterCellValue(a,l)
  w = apply(sum,v)
  o = [ sum(vi) for vi in v ] 
  test_iter_cell_value(w,o)
end

for a in ax
  v = TestIndexCellValue(a,l)
  w = apply(-,v)
  o = [ -vi for vi in v ] 
  test_index_cell_value(w,o)
end

for (a,b) in zip(ax,bx)
  u = TestIndexCellValue(a,l)
  v = TestIndexCellValue(b,l)
  w = apply(-,u,v)
  o = [ ui-vi for (ui,vi) in zip(u,v) ] 
  test_index_cell_value(w,o)
end

for a in [[1,2,3], SVector(1.3,2.0), VectorValue(1.3,2.0)]
  v = TestIndexCellValue(a,l)
  w = apply(sum,v)
  o = [ sum(vi) for vi in v ] 
  test_index_cell_value(w,o)
end

end # module
