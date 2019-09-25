module VtkioTests

using Test
using Gridap

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

  grid = CartesianGrid(partition=(3,3))
  trian = Triangulation(grid)
  quad = CellQuadrature(trian,degree=2)

  phi = CellGeomap(trian)

  q = coordinates(quad)

  x = evaluate(phi,q)

  ufun(x) = 2*x[1] + x[2]
  u = CellField(trian,ufun)

  vfun(x) = VectorValue(x[1],1.0)
  v = CellField(trian,vfun)

  tfun(x) = TensorValue(x[1],1.0,x[2],0.1)
  t = CellField(trian,tfun)

  pdata =[
    "u"=>evaluate(u,q),"v"=>evaluate(v,q),
    "t"=>evaluate(t,q),"t2"=>apply(tfun,x,broadcast=true)]

  writevtk(x,f)
  writevtk(x,f,pointdata=pdata)

  rm(d,recursive=true)

end

@testset "WritevtkForCellPoint" begin

  d = mktempdir()
  f = joinpath(d,"x")

  grid = CartesianGrid(partition=(3,3))
  trian = Triangulation(grid)

  xe = CellPoints(trian)

  x = cellmean(xe)

  tfun(x) = TensorValue(x[1],1.0,x[2],0.1)

  t = apply(tfun,x)

  writevtk(x,f,pointdata=["t"=>t])

  rm(d,recursive=true)

end

@testset "WritePolytope" begin

  d = mktempdir()
  f = joinpath(d,"polytope")

  polytope = Polytope(1,1,1)

  writevtk(polytope,f)

  rm(d,recursive=true)

end

@testset "CartesianDiscreteModel" begin

  d = mktempdir()
  f = joinpath(d,"model")
  #f = "model"

  model = CartesianDiscreteModel(
    domain=(0,1,0,1,0,1),
    partition=(2,2,2))

  writevtk(model,f)

  rm(d,recursive=true)

end

@testset "WritevtkForIntegrationMesh" begin

grid = CartesianGrid(partition=(3,3))
trian = Triangulation(grid)

ufun(x) = 3*x[2]*x[1]
u = CellField(trian,ufun)

vfun(x) = VectorValue(3*x[2]*x[1],2*x[2])
v = CellField(trian,vfun)

d = mktempdir()
f = joinpath(d,"trian")

writevtk(trian,f)
writevtk(trian,f,nref=2,)
writevtk(trian,f,nref=2,celldata=["r"=>rand(9)])
writevtk(trian,f,nref=2,cellfields=["u"=>u,"v"=>v])
writevtk(trian,f,nref=2,celldata=["r"=>rand(9)],cellfields=["u"=>u,"v"=>v])

grid = CartesianGrid(partition=(3,3))
grid = simplexify(grid)
trian = Triangulation(grid)
writevtk(trian,f)

rm(d,recursive=true)

end

end # module VtkioTests
