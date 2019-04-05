module CellIntegrationTests

using Test
using Numa.CellValues
using Numa.CellFunctions
using Numa.CellQuadratures
using Numa.CellIntegration
using Numa.Quadratures

include("CellIntegrationTestsMocks.jl")

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

@testset "Integrate" begin

  basis = cellbasis(imesh)

  refquad = TensorProductQuadrature(orders=(2,2))

  quad = ConstantCellQuadrature(refquad,ncells(imesh))

  mmat = integrate(inner(basis,basis),imesh,quad)

  @test isa(mmat,CellArray{Float64,2})

  ufun(x::Point{2}) = 1.0

  cellvol = integrate(ufun,imesh,quad)

  @test isa(cellvol,CellValue{Float64})

  for vi in cellvol
    @assert vi ≈ (1.0/3)^2
  end

end

end # module CellIntegrationTests
