module CartesianTests

using Test
using Numa.Geometry
using Numa.Geometry: CartesianGridCoords
using Numa.Geometry: CartesianGridConnectivity

@testset "Cartesian" begin

  domain = ((0.0,1.0),(-1.0,2.0),(2.0,3.0))
  npoints = (3,2,4)
  x = CartesianGridCoords(domain,npoints)

  for xi in x
    @show xi
  end

  ncells = (2,3)
  t = CartesianGridConnectivity(ncells)

  for ti in t
    @show ti
  end

  @show length(t)
  @show size(t)

end

end # module CartesianTests
