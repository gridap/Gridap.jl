module CellIntegrationTests

using Test
using Numa
using Numa.FieldValues
using Numa.CellValues
using Numa.CellFunctions
using Numa.CellQuadratures
using Numa.CellIntegration
using Numa.Quadratures
using Numa.Geometry
using Numa.Geometry.Cartesian

grid = CartesianGrid(partition=(3,3))
trian = triangulation(grid)

@testset "Geomap" begin

  phi = geomap(trian)

  @test isa(phi,CellGeomap{2,2})

end

@testset "Integrate" begin

  basis = cellbasis(trian)

  refquad = TensorProductQuadrature(orders=(2,2))

  quad = ConstantCellQuadrature(refquad,ncells(trian))

  mmat = integrate(inner(basis,basis),trian,quad)

  @test isa(mmat,CellArray{Float64,2})

  ufun(x::Point{2}) = 1.0

  cellvol = integrate(ufun,trian,quad)

  @test isa(cellvol,CellValue{Float64})

  for vi in cellvol
    @assert vi ≈ (2.0/3)^2
  end

end

@testset "Cellfield" begin

  ufun(x::Point{2}) = 2*x[1]+x[2]

  u = cellfield(trian,ufun)

  νfun(x::Point{2},u::Float64) = TensorValue(x[1], u*x[2], 0.0, u)

  ν(u) = cellfield(trian,νfun,u)

  @test isa(u,CellField{2,Float64})
  @test isa(ν(u),CellField{2,TensorValue{2,4}})

end

end # module CellIntegrationTests
