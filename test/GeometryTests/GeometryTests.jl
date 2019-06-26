module GeometryTests

using Test
using Gridap

@testset "CartesianGrid" begin

  grid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

  @test celldim(grid) == 2
  @test pointdim(grid) == 2
  @test ncells(grid) == 3*4
  @test npoints(grid) == (3+1)*(4+1)

  x = points(grid)

  @test isa(x,CellValue{Point{2,Float64}})

  @test length(x) == 20
  @test x[13] == Point(0.0,1.25)
  @test x[4] == Point(1.0,-1.0)
  @test x[1] == Point(0.0,-1.0)
  @test x[end] == Point(1.0,2.0)

  t = cells(grid)

  @test isa(t,CellArray{Int,1})

  @test length(t) == 12

  @test t[1] == [1,2,5,6]
  @test t[11] == [14,15,18,19]

  c = celltypes(grid)
  @test isa(c,ConstantCellValue{NTuple{2,Int}})
  @test c.value == (HEX_AXIS,HEX_AXIS)
  @test length(c) == 12

  o = cellorders(grid)
  @test isa(o,ConstantCellValue{Int})
  @test o.value == 1
  @test length(o) == 12

  #graph = gridgraph(grid)

  #@test isa(graph,GridGraph)

  #@test isa(celltovefs(graph), IndexCellArray{Int,1})

  #@test isa(veftocells(graph), IndexCellArray{Int,1})

  trian = Triangulation(grid)

  xe = CellPoints(trian)
  @test isa(xe,CellPoints{2,Float64})

  cb = CellBasis(trian)
  @test isa(cb,CellBasis{2,Float64})

end

@testset "FlexibleUnstructuredGrid" begin

  cgrid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

  grid = FlexibleUnstructuredGrid(cgrid)

  c = celltypes(grid)

  @test isa(c,ConstantCellValue{NTuple{2,Int}})
  @test c.value == (HEX_AXIS,HEX_AXIS)
  @test length(c) == 12

  o = cellorders(grid)
  @test isa(o,ConstantCellValue{Int})
  @test o.value == 1
  @test length(o) == 12

  #d = mktempdir()
  #f = joinpath(d,"grid")

  #writevtk(grid,f)

  #rm(d,recursive=true)

end

@testset "UnstructuredGrid" begin

  cgrid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))

  grid = UnstructuredGrid(cgrid)

  c = celltypes(grid)
  @test isa(c,ConstantCellValue{NTuple{2,Int}})
  @test c.value == (HEX_AXIS,HEX_AXIS)
  @test length(c) == 12

  o = cellorders(grid)
  @test isa(o,ConstantCellValue{Int})
  @test o.value == 1
  @test length(o) == 12

  #d = mktempdir()
  #f = joinpath(d,"grid")

  #writevtk(grid,f)

  #rm(d,recursive=true)

end

@testset "FaceLabels" begin

  vertex_to_geolabel = [1,1,2,2,2,1,1,3,3]
  edge_to_geolabel = [4,4,5,5,5,5,6,6,4]
  physlabel_1 = [1,3,4]
  physlabel_2 = [5,3,6,2]
  tag_to_name = ["label1","label2"]

  labels = FaceLabels(
    [vertex_to_geolabel, edge_to_geolabel],
    [physlabel_1, physlabel_2],
    tag_to_name)

  @test isa(labels,FaceLabels)
  @test labels_on_dim(labels,0) == vertex_to_geolabel
  @test labels_on_dim(labels,1) == edge_to_geolabel
  @test labels_on_tag(labels,1) == physlabel_1
  @test labels_on_tag(labels,2) == physlabel_2
  @test tag_from_name(labels,"label1")==1

end

@testset "CartesianDiscreteModel" begin

  model = CartesianDiscreteModel(
    domain=(0.0,1.0,-1.0,2.0,0.0,1.0),
    partition=(3,4,2))

  grid3 = Grid(model,3)

  @test isa(grid3,CartesianGrid{3})

  gridgraph = FullGridGraph(model)

  @test isa(gridgraph, FullGridGraph)

  grid2 = Grid(model,2)
  grid1 = Grid(model,1)
  grid0 = Grid(model,0)

  labels = FaceLabels(model)

  @test isa(labels,FaceLabels)
  @test length(labels_on_dim(labels,3)) == ncells(grid3)
  @test length(labels_on_dim(labels,2)) == ncells(grid2)
  @test length(labels_on_dim(labels,1)) == ncells(grid1)
  @test length(labels_on_dim(labels,0)) == ncells(grid0)

  @test ntags(labels) == 28
  @test name_from_tag(labels,28) == "boundary"

end

end # module GeometryTests
