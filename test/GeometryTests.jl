module GeometryTests

using Test
using Numa
using Numa.FieldValues
using Numa.CellValues
using Numa.CellMaps
using Numa.Geometry
using Numa.Geometry.Cartesian
using Numa.Geometry.Unstructured
using Numa.Geometry.Wrappers
using Numa.Polytopes
using Numa.Vtkio

@testset "CartesianGrid" begin

  grid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

  x = points(grid)

  @test isa(x,CellValue{Point{2}})

  @test length(x) == 20
  @test x[13] == [0.0,1.25]
  @test x[4] == [1.0,-1.0]
  @test x[1] == [0.0,-1.0]
  @test x[end] == [1.0,2.0]

  t = cells(grid)

  @test isa(t,CellArray{Int,1})

  @test length(t) == 12

  @test t[1] == [1,2,5,6]
  @test t[11] == [14,15,18,19]

  c = celltypes(grid)
  @test isa(c,ConstantCellValue{NTuple{2,Int}})
  @test celldata(c) == (HEX_AXIS,HEX_AXIS)
  @test length(c) == 12

  o = cellorders(grid)
  @test isa(o,ConstantCellValue{Int})
  @test celldata(o) == 1
  @test length(o) == 12

  graph = gridgraph(grid)

  @test isa(graph,GridGraph)

  @test isa(celltovefs(graph), IndexCellArray{Int,1})

  @test isa(veftocells(graph), IndexCellArray{Int,1})

  trian = triangulation(grid)

  xe = cellcoordinates(trian)
  @test isa(xe,CellPoints{2})

  cb = cellbasis(trian)
  @test isa(cb,CellBasis{2,ScalarValue})

end

@testset "FlexibleUnstructuredGrid" begin

  cgrid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

  grid = FlexibleUnstructuredGrid(cgrid)

  c = celltypes(grid)

  @test isa(c,ConstantCellValue{NTuple{2,Int}})
  @test celldata(c) == (HEX_AXIS,HEX_AXIS)
  @test length(c) == 12

  o = cellorders(grid)
  @test isa(o,ConstantCellValue{Int})
  @test celldata(o) == 1
  @test length(o) == 12

  d = mktempdir()
  f = joinpath(d,"grid")

  writevtk(grid,f)

  rm(d,recursive=true)

end

@testset "UnstructuredGrid" begin

  cgrid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

  grid = UnstructuredGrid(cgrid)

  c = celltypes(grid)
  @test isa(c,ConstantCellValue{NTuple{2,Int}})
  @test celldata(c) == (HEX_AXIS,HEX_AXIS)
  @test length(c) == 12

  o = cellorders(grid)
  @test isa(o,ConstantCellValue{Int})
  @test celldata(o) == 1
  @test length(o) == 12

  d = mktempdir()
  f = joinpath(d,"grid")

  writevtk(grid,f)

  rm(d,recursive=true)

end

@testset "NFacesLabels" begin

  vertex_to_geolabel = [1,1,2,2,2,1,1,3,3]
  edge_to_geolabel = [4,4,5,5,5,5,6,6,4]
  physlabel_1 = [1,3,4]
  physlabel_2 = [5,3,6,2]
  
  nfacelabels = NFacesLabels(
    (vertex_to_geolabel, edge_to_geolabel),
    [physlabel_1, physlabel_2])
  
  @test isa(nfacelabels,NFacesLabels{1})
  @test nfacegeolabel(nfacelabels,0) == vertex_to_geolabel
  @test nfacegeolabel(nfacelabels,1) == edge_to_geolabel
  @test geolabels(nfacelabels,1) == physlabel_1
  @test geolabels(nfacelabels,2) == physlabel_2

end

@testset "GridGraph" begin

  # Dummy data. Not related with an actual grid
  cell_to_vertices = CellArrayFromArrayOfArrays([[1,2,4,2], [2,3,1,3]])
  cell_to_edges = CellArrayFromArrayOfArrays([[1,2,3,4], [4,1,2,3]])
  vertex_to_cells = CellArrayFromArrayOfArrays([[1,2],[2],[1],[1]])
  edge_to_cells = CellArrayFromArrayOfArrays([[1],[1,2],[2,1],[2]])

  graph = NewGridGraph(
    (cell_to_vertices, cell_to_edges),
    (vertex_to_cells, edge_to_cells))

  @test isa(graph,NewGridGraph{1})
  @test cellvefs(graph,0) == cell_to_vertices
  @test cellvefs(graph,1) == cell_to_edges
  @test vefcells(graph,0) == vertex_to_cells
  @test vefcells(graph,1) == edge_to_cells

end

@testset "RefCellFromPolytope" begin

  polytope = Polytope((1,1,1))

  refcell = RefCell(polytope)

end

end # module GeometryTests
