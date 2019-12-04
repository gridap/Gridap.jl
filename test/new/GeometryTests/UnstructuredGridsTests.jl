module UnstructuredGridsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry: ConformingTrianMock

# Unstructured grid from raw data

trian = ConformingTrianMock()

node_coordinates = get_node_coordinates(trian)
cell_nodes = get_cell_nodes(trian)
reffes = get_reffes(trian)
cell_types = get_cell_type(trian)

grid = UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types)
test_conforming_triangulation(grid)

q1i = Point(0.5,0.5)
np1 = 4
q1 = fill(q1i,np1)

q2i = Point(0.25,0.25)
np2 = 3
q2 = fill(q2i,np2)

q = CompressedArray([q1,q2],get_cell_type(grid))

cell_map = get_cell_map(grid)
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


# UnstructuredGrid from ConformingTriangulation

grid = UnstructuredGrid(trian)
test_conforming_triangulation(grid)
@test grid === UnstructuredGrid(grid)

# get low dim grid

grid1 = UnstructuredGrid(ReferenceFE{1},trian)
@test num_cell_dims(grid1) == 1

# UnstructuredGrid from NodalReferenceFE

quad8 = LagrangianRefFE(Float64,QUAD,2)
grid = UnstructuredGrid(quad8)
@test num_nodes(grid) == num_nodes(quad8)

# UnstructuredGrid from Polytope

grid = UnstructuredGrid(ReferenceFE{2},WEDGE)
@test num_cells(grid) == 5
@test num_cell_dims(grid) == 2
@test num_point_dims(grid) == 3

grid = UnstructuredGrid(ReferenceFE{3},WEDGE)
@test num_cells(grid) == 1
@test num_cell_dims(grid) == 3
@test num_point_dims(grid) == 3

end # module
