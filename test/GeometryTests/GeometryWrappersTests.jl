module GeometryWrappersTests

using Gridap
using UnstructuredGrids: RefCell, UGrid
using Test
using UnstructuredGrids: generate_dual_connections
using UnstructuredGrids: find_cell_to_faces

# Polytope to UnstructuredGrid

t = TET_AXIS
polytope = Polytope((t,t,t))

grid = Grid(polytope,1)
grid = Grid(polytope,2)
grid = Grid(polytope,3)
#writevtk(grid,"grid")

t = HEX_AXIS
polytope = Polytope((t,t,t))

grid = Grid(polytope,1)
grid = Grid(polytope,2)
grid = Grid(polytope,3)
#writevtk(grid,"grid")

trian = Triangulation(grid)
quad = CellQuadrature(trian,order=2)

q = coordinates(quad)
phi = CellGeomap(trian)
x = evaluate(phi,q)

#writevtk(x,"x")

# Polytope to RefCell

t = HEX_AXIS
polytope = Polytope((t,t,t))

refcell = RefCell(polytope)

# UnstructuredGrid to UGrid

c2v = [[1,2,4,5],[2,3,5,6]]
v2x = Point{2,Float64}[(0.,0.),(0.,1.),(0.,2.),(1.,0.),(1.,1.),(1.,2.)]
l = length(c2v)
order = 1
t = (HEX_AXIS,HEX_AXIS)
c2t = ConstantCellValue(t,l)
c2o = ConstantCellValue(order,l)
c2v_data, c2v_ptrs = generate_data_and_ptrs(c2v)

grid = UnstructuredGrid(v2x,c2v_data,c2v_ptrs,c2t,c2o)
ugrid = UGrid(grid)

graph = FullGridGraph(grid)
test_full_grid_graph(graph,2)

t = TET_AXIS
polytope = Polytope((t,t,t))

grid = Grid(polytope,2)

# Others

model = CartesianDiscreteModel(partition=(4,4,3))
graph = FullGridGraph(model)
edge_to_vertices = connections(graph,1,0)
vertex_to_edges = connections(graph,0,1)
r = generate_dual_connections(edge_to_vertices)
@test r == vertex_to_edges

grid2 = Grid(model,2)
vertex_to_edges = connections(graph,0,1)
face_to_edges = connections(graph,2,1)
r = find_cell_to_faces(grid2,vertex_to_edges,1)
@test r == face_to_edges


end # module

