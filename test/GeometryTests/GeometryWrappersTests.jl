module AAA

using Gridap
using Gridap.Helpers
using UnstructuredGrids.Kernels: generate_data_and_ptrs
using Gridap.GeometryWrappers: _faces

import Gridap: Grid

function Grid(polytope::Polytope{D},dim::Int) where D
  @assert dim < D

  orders = fill(1,D)
  na = NodesArray(polytope,orders)

  dim_to_jface_to_vertices, dim_to_jface_to_code = _faces(polytope)
  jface_to_vertices = dim_to_jface_to_vertices[dim+1]
  jface_to_code = dim_to_jface_to_code[dim+1]

  njfaces = length(jface_to_code)
  @assert njfaces > 0
  code1 = jface_to_code[1]
  @notimplementedif any([ code1 != code for code in jface_to_code ])

  points = na.coordinates
  cells_data, cells_ptrs = generate_data_and_ptrs(jface_to_vertices)
  code = (code1...,)
  ctypes = ConstantCellValue(code,njfaces)
  order = 1
  orders = ConstantCellValue(order,njfaces)

  UnstructuredGrid(points,cells_data,cells_ptrs,ctypes,orders)

end

end # module


module GeometryWrappersTests

using Gridap
using UnstructuredGrids: RefCell, UGrid
using Test

# Polytope to UnstructuredGrid

t = HEX_AXIS
polytope = Polytope((t,t,t))

grid = Grid(polytope,2)
writevtk(grid,"grid")

#trian = Triangulation(grid)
#quad = CellQuadrature(trian,order=5)
#
#q = coordinates(quad)
#phi = CellGeomap(trian)
#x = evaluate(phi,q)

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

