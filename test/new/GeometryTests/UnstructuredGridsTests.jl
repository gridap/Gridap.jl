module UnstructuredGridsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry: ConformingTrianMock

node_coordinates = Point{2,Float64}[(0,0),(1,0),(2,0),(1,1),(2,1),(0,2),(2,2)]

cell_nodes = Table([[1,2,6,4],[2,3,4,5],[4,5,7],[4,7,6]])

order = 1
tri3 = LagrangianRefFE(Float64,TRI,order)
quad4 = LagrangianRefFE(Float64,QUAD,order)
reffes = [tri3, quad4]

cell_types = [2,2,1,1]

trian = UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types)

test_conforming_triangulation(trian)

q1i = Point(0.25,0.25)
np1 = 3
q1 = fill(q1i,np1)
q2i = Point(0.5,0.5)
np2 = 4
q2 = fill(q2i,np2)
q = CompressedArray([q1,q2],get_cell_types(trian))

cell_map = get_cell_map(trian)
x = evaluate(cell_map,q)

x1i = Point(0.5, 0.75)
x2i = Point(1.5, 0.5)
x3i = Point(1.5, 1.25)
x4i = Point(1.0, 1.5)
x1 = fill(x1i,np2)
x2 = fill(x2i,np2)
x3 = fill(x3i,np1)
x4 = fill(x4i,np1)
x = [x1,x2,x3,x4]
test_array_of_fields(cell_map,q,x)

trian = ConformingTrianMock()

grid = UnstructuredGrid(trian)
test_conforming_triangulation(trian)

@test grid === UnstructuredGrid(grid)

end # module
