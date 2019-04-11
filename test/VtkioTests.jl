module VtkioTests

using Test
using Numa
using Numa.Quadratures
using Numa.CellQuadratures
using Numa.CellFunctions
using Numa.CellValues
using Numa.Geometry
using Numa.Polytopes
using Numa.Vtkio

include("CellIntegrationTestsMocks.jl")

@testset "VTKioGrid" begin

  d = mktempdir()
  f = joinpath(d,"grid")

  grid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0,0.0,1.0),partition=(10,10,10))

  cd1 = rand(length(cells(grid)))
  cd2 = 1:length(cells(grid))

  pd1 = rand(length(points(grid)))
  pd2 = 1:length(points(grid))
  pd3 = [ VectorValue(1.0,2.0) for i in 1:length(points(grid)) ]

  cdat = ["cd1"=>cd1,"cd2"=>cd2]
  pdat = ["pd1"=>pd1,"pd2"=>pd2,"pd3"=>pd3]

  writevtk(grid,f)
  writevtk(grid,f,celldata=cdat)
  writevtk(grid,f,pointdata=pdat)
  writevtk(grid,f,celldata=cdat,pointdata=pdat)

  rm(d,recursive=true)

end

@testset "WritevtkForCellPoints" begin

  d = mktempdir()
  f = joinpath(d,"x")

  imesh = DummyIntegrationMesh2D(partition=(3,3))
  refquad = TensorProductQuadrature(orders=(2,2))
  quad = ConstantCellQuadrature(refquad,ncells(imesh))

  phi = geomap(imesh)

  q = coordinates(quad)

  x = evaluate(phi,q)

  ufun(x) = 2*x[1] + x[2]
  u = cellfield(imesh,ufun)

  vfun(x) = VectorValue(x[1],1.0)
  v = cellfield(imesh,vfun)

  tfun(x) = TensorValue(x[1],1.0,x[2],0.1)
  t = cellfield(imesh,tfun)

  pdata =["u"=>evaluate(u,q),"v"=>evaluate(v,q),"t"=>evaluate(t,q),"t2"=>apply(tfun,x)]

  writevtk(x,f)
  writevtk(x,f,pointdata=pdata)

  rm(d,recursive=true)

end

@testset "WritevtkForCellPoint" begin

  d = mktempdir()
  f = joinpath(d,"x")

  imesh = DummyIntegrationMesh2D(partition=(3,3))

  xe = cellcoordinates(imesh)

  x = cellmean(xe)

  tfun(x) = TensorValue(x[1],1.0,x[2],0.1)

  t = apply(tfun,x)

  writevtk(x,f,pointdata=["t"=>t])

  rm(d,recursive=true)

end

@testset "VisualizationGrid" begin

  imesh = DummyIntegrationMesh2D(partition=(3,3))

  vg = visgrid(imesh,nref=2)

  ufun(x) = 3*x[2]*x[1]
  u = cellfield(imesh,ufun)

  vfun(x) = VectorValue(3*x[2]*x[1],2*x[2])
  v = cellfield(imesh,vfun)

  d = mktempdir()
  f = joinpath(d,"grid")

  writevtk(vg,f,)
  writevtk(vg,f,celldata=["r"=>rand(9)])
  writevtk(vg,f,cellfields=["u"=>u,"v"=>v])
  writevtk(vg,f,celldata=["r"=>rand(9)],cellfields=["u"=>u,"v"=>v])

  rm(d,recursive=true)

end

end # module VtkioTests
