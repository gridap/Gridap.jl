
@testset "CachedArray" begin

  using Numa.CellValues: CachedArray, setsize!

  x = rand(10,5)

  y = CachedArray(x)

  @test y == x

  z = y + x

  a = y.array

  setsize!(y,(5,3))

  @test y.array === a

  @test y == x[1:5,1:3]

  setsize!(y,(15,30))

  @test ! ( y.array === a )

  @test size(y) == size(y.array)


end
