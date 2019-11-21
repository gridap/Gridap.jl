module HelpersTests

using Test
using Gridap.Helpers

for D = 0:5
  @test tfill(2,Val(D)) == tuple(fill(2,D)...)
end

end # module

