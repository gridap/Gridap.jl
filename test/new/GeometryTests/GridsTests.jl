module GridsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs

using Gridap.Geometry: GridMock

trian = GridMock()
test_grid(trian)

q1i = Point(0.5,0.5)
np1 = 4
q1 = fill(q1i,np1)

q2i = Point(0.25,0.25)
np2 = 3
q2 = fill(q2i,np2)

q = CompressedArray([q1,q2],get_cell_type(trian))

cell_map = get_cell_map(trian)
x = evaluate(cell_map,q)

x1i = Point(0.5, 0.5)
x2i = Point(1.25, 0.25)
x3i = Point(1.75, 0.5)
x4i = Point(0.5, 1.5)
x5i = Point(1.5, 1.5)
x1 = fill(x1i,np1)
x2 = fill(x2i,np2)
x3 = fill(x3i,np2)
x4 = fill(x4i,np1)
x5 = fill(x5i,np1)
x = [x1,x2,x3,x4,x5]

test_array_of_fields(cell_map,q,x)

# from LagrangianRefFE

quad8 = LagrangianRefFE(Float64,QUAD,2)
grid = Grid(quad8)
@test num_nodes(grid) == num_nodes(quad8)

# from Polytope

grid = Grid(ReferenceFE{2},WEDGE)
@test num_cells(grid) == 5
@test num_cell_dims(grid) == 2
@test num_point_dims(grid) == 3

grid = Grid(ReferenceFE{3},WEDGE)
@test num_cells(grid) == 1
@test num_cell_dims(grid) == 3
@test num_point_dims(grid) == 3

## get low dim grid
#
#grid = Grid(ReferenceFE{1},trian)

end # module
