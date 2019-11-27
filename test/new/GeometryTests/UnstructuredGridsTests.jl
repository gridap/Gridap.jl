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
cell_types = get_cell_types(trian)

grid = UnstructuredGrid(node_coordinates,cell_nodes,reffes,cell_types)
test_conforming_triangulation(grid)

q1i = Point(0.5,0.5)
np1 = 4
q1 = fill(q1i,np1)

q2i = Point(0.25,0.25)
np2 = 3
q2 = fill(q2i,np2)

q = CompressedArray([q1,q2],get_cell_types(grid))

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


# UnstructuredGrid from NodalReferenceFE

quad8 = LagrangianRefFE(Float64,QUAD,2)
grid = UnstructuredGrid(quad8)
@test num_nodes(grid) == num_nodes(quad8)


import Gridap.Geometry: UnstructuredGrid

function UnstructuredGrid(
  ::Type{ReferenceFE{d}},
  grid::UnstructuredGrid,
  cell_to_nodes::Table,
  cell_to_faces::Table) where d

  ctype_to_reffe = get_reffes(grid)
  cell_to_cell_type = get_cell_types(grid)

  t = _generate_ftype_to_refface(Val{d}(),ctype_to_reffe)
  ftype_to_refface, ctype_to_lface_to_ftype = t

  face_to_ftype = generate_face_to_face_type(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype)

  ctype_to_lface_to_lnodes = map(
    (r) -> get_face_nodeids(r)[d+1] ,ctype_to_reffe)

  face_to_nodes = generate_face_to_vertices(
    cell_to_nodes,
    cell_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lnodes)

  node_to_coord = get_node_coordinates(grid)

  UnstructuredGrid(node_to_coord,face_to_nodes,ftype_to_refface,face_to_ftype)

end

## UnstructuredGridGraph
#
#using Gridap.Helpers
#using Gridap.Arrays
#include("../../../src/new/Geometry/GridOperations.jl")
#
#struct UnstructuredGridGraph <: GridGraph
#  d_to_cell_to_dfaces::Vector{Table{Int,Int32}}
#  d_to_dface_to_cells::Vector{Table{Int,Int32}}
#  vertex_to_node::Vector{Int}
#end
#
#function UnstructuredGridGraph(trian::ConformingTriangulation)
#  grid = UnstructuredGrid(trian)
#  UnstructuredGridGraph(grid)
#end
#
#"""
#"""
#function UnstructuredGridGraph(grid::UnstructuredGrid)
#
#  primal, dual, vertex_to_node = _generate_grid_graph(grid)
#  UnstructuredGridGraph(primal, dual, vertex_to_node)
#end
#
#function _generate_grid_graph(grid)
#
#  D = num_cell_dims(grid)
#  cell_to_vertices, vertex_to_node = generate_cell_to_vertices(grid)
#  vertex_to_cells = generate_cells_around(cell_to_vertices)
#
#  T = typeof(cell_to_vertices)
#  primal = Vector{T}(undef,D+1)
#  dual = Vector{T}(undef,D)
#
#  d = 0
#  primal[d+1] = cell_to_vertices
#  dual[d+1] = vertex_to_cells
#
#  for d in 1:(D-1)
#    cell_to_dfaces = generate_cell_to_faces( d, grid, cell_to_vertices, vertex_to_cells)
#    dfaces_to_cells = generate_cells_around(cell_to_dfaces)
#    primal[d+1] = cell_to_dfaces
#    dual[d+1] = dfaces_to_cells
#  end
#
#  ncells = length(cell_to_vertices)
#  primal[D+1] = identity_table(Int,Int32,ncells)
#
#
#  (primal, dual, vertex_to_node)
#end
#
#graph = UnstructuredGridGraph(grid)
##display(graph.d_to_cell_to_dfaces)
##display(graph.d_to_dface_to_cells)
##display(graph.vertex_to_node)
#
##test_grid_graph(graph)

end # module
