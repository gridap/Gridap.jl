
@testset "CachedArray" begin

  using Numa.CellValues: CachedArray

  x = rand(10,5)

  y = CachedArray(x)

  @test y == x

  z = y + x

end
