module HelpersTests

using Test
using Gridap.Helpers

@testset "PtrsUtils" begin

  a = [9,2,1,2,4,7,4]
  b = [1,9,2,1,2,4,7]

  rewind_ptrs!(a)
  @test a == b

  a = [3,2,4,2]
  b = [1,3,7,9]

  length_to_ptrs!(a)
  @test a == b

end

end # module HelpersTests
