module CellValuesMocksTests

using Test
using Gridap
using StaticArrays
using TensorValues

using ..CellValuesMocks

l = 10

for a in [1.2, SVector(1.3,2.0), VectorValue(1.3,2.0)]

  v = TestIterCellValue(a,l)
  
  w = [a for i in 1:l]
  test_iter_cell_value(v,w)
  
  v = TestIndexCellValue(a,l)
  test_index_cell_value(v,w)

end

a = 1.2
v = TestIndexCellValue(a,l)

s = string(v)
z = """
1 -> 1.2
2 -> 1.2
3 -> 1.2
4 -> 1.2
5 -> 1.2
6 -> 1.2
7 -> 1.2
8 -> 1.2
9 -> 1.2
10 -> 1.2
"""
@test s == z


end # module
