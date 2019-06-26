module GeometryWrappersTests

using Gridap
using UnstructuredGrids: RefCell, UGrid
using Test

# Polytope to UnstructuredGrid

t = HEX_AXIS
polytope = Polytope((t,t,t))

grid = Grid(polytope,2)
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

# CartesianGrid to UnstructuredGrid

cgrid = CartesianGrid(
  domain=(0.0,1.0,-1.0,2.0,0.0,1.0),
  partition=(3,4,2))

# UnstructuredGrid to UGrid

grid = UnstructuredGrid(cgrid)

ugrid = UGrid(grid)

graph = FullGridGraph(grid)
@test isa(connections(graph,3,0), CellArray)
@test isa(connections(graph,3,1), CellArray)
@test isa(connections(graph,3,2), CellArray)
@test isa(connections(graph,2,0), CellArray)
@test isa(connections(graph,2,1), CellArray)
@test isa(connections(graph,2,3), CellArray)
@test isa(connections(graph,1,0), CellArray)
@test isa(connections(graph,1,2), CellArray)
@test isa(connections(graph,1,3), CellArray)
@test isa(connections(graph,0,0), CellArray)
@test isa(connections(graph,0,1), CellArray)
@test isa(connections(graph,0,2), CellArray)
@test isa(connections(graph,0,3), CellArray)

end # module

