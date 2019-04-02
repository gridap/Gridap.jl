module IntegrationMeshesTests

using Test
using Numa.CellFunctions
using Numa.IntegrationMeshes

include("IntegrationMeshesTestsMocks.jl")

imesh = DummyIntegrationMesh2D(partition=(3,3))

@testset "Mocks" begin

  @test isa(imesh,IntegrationMesh)
  
  coords = cellcoordinates(imesh)
  
  basis = cellbasis(imesh)
  
  phi = geomap(imesh)
  
  @test isa(coords,CellPoints{2})
  
  @test isa(basis,CellBasis{2,Float64})
  
  @test isa(phi,CellField{2,Point{2}})

end

@testset "Geomap" begin

  phi = geomap(imesh)

  @test isa(phi,CellGeomap{2,2})

end

end # module IntegrationMeshesTests
