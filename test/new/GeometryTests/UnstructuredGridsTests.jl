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

quad8 = LagrangianRefFE(Float64,QUAD,2)

grid = UnstructuredGrid(quad8)
@test num_nodes(grid) == num_nodes(quad8)

# UnstructuredGridGraph

using Gridap.Helpers
using Gridap.Arrays
include("../../../src/new/Geometry/GridOperations.jl")

struct UnstructuredGridGraph <: GridGraph
  d_to_cell_to_dfaces::Vector{Table{Int,Int32}}
  d_to_dface_to_cells::Vector{Table{Int,Int32}}
  vertex_to_node::Vector{Int}
end

function UnstructuredGridGraph(trian::ConformingTriangulation)
  grid = UnstructuredGrid(trian)
  UnstructuredGridGraph(grid)
end

"""
"""
function UnstructuredGridGraph(grid::UnstructuredGrid)

  primal, dual, vertex_to_node = _generate_grid_graph(grid)
  UnstructuredGridGraph(primal, dual, vertex_to_node)
end

function _generate_grid_graph(grid)

  D = num_cell_dims(grid)
  cell_to_vertices, vertex_to_node = generate_cell_to_vertices(grid)
  vertex_to_cells = generate_cells_around(cell_to_vertices)

  T = typeof(cell_to_vertices)
  primal = Vector{T}(undef,D+1)
  dual = Vector{T}(undef,D)

  d = 0
  primal[d+1] = cell_to_vertices
  dual[d+1] = vertex_to_cells

  for d in 1:(D-1)
    cell_to_dfaces = generate_cell_to_faces( d, grid, cell_to_vertices, vertex_to_cells)
    dfaces_to_cells = generate_cells_around(cell_to_dfaces)
    primal[d+1] = cell_to_dfaces
    dual[d+1] = dfaces_to_cells
  end

  ncells = length(cell_to_vertices)
  primal[D+1] = identity_table(Int,Int32,ncells)


  (primal, dual, vertex_to_node)
end

graph = UnstructuredGridGraph(grid)
#display(graph.d_to_cell_to_dfaces)
#display(graph.d_to_dface_to_cells)
#display(graph.vertex_to_node)

#test_grid_graph(graph)

end # module
