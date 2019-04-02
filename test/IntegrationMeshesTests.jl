module IntegrationMeshesTests

using Test
using Numa.CellArrays
using Numa.CellFunctions
using Numa.CellQuadratures
using Numa.IntegrationMeshes
using Numa.Quadratures

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

@testset "Integrate" begin

  basis = cellbasis(imesh)

  refquad = TensorProductQuadrature(orders=(2,2))

  quad = ConstantCellQuadrature(refquad,ncells(imesh))

  mmat = integrate(inner(basis,basis),imesh,quad)

  @test isa(mmat,CellArray{Float64,2})

  @eval begin
    ufun(x::Point{2}) = 1
    ufun(::Type{Point{2}}) = Float64
  end

  cellvol = integrate(ufun,imesh,quad)

  # @fverdugo TODO: it should return an array scalar
  @test isa(cellvol,CellArray{Float64,0})

  for vi in cellvol
    @assert vi[1] ≈ (1.0/3)^2
  end

end

end # module IntegrationMeshesTests
