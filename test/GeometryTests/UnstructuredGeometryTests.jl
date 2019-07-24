module UnstructuredGeometryTests

using Test
using Gridap

c2v = [[1,2,4,5],[2,3,5,6]]
v2x = Point{2,Float64}[(0.,0.),(0.,1.),(0.,2.),(1.,0.),(1.,1.),(1.,2.)]

l = length(c2v)
order = 1
t = (HEX_AXIS,HEX_AXIS)
c2t = ConstantCellValue(t,l)
c2o = ConstantCellValue(order,l)

c2v_data, c2v_ptrs = generate_data_and_ptrs(c2v)

grid = UnstructuredGrid(v2x,c2v_data,c2v_ptrs,c2t,c2o)
@test isa(grid,UnstructuredGrid)
test_grid(grid,6,2)

grid = FlexibleUnstructuredGrid(v2x,c2v,c2t,c2o)
@test isa(grid,FlexibleUnstructuredGrid)
test_grid(grid,6,2)

end # module
