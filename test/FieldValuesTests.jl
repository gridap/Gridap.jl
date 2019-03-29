module FieldValuesTests

using Test
using Numa.FieldValues

@testset "FieldValuesConstructors" begin

  vv = VectorValue(1.0,2.0,2.3)

  @test isa(vv,VectorValue{3})

  mvv = MVectorValue(1.0,2.0,2.3)

  @test isa(mvv,MVectorValue{3})

  tv = TensorValue(1.0,2.0,2.3,2.3)

  @test isa(tv,TensorValue{2,4})
  @test tv[2,1] == 2.0

  mtv = MTensorValue(1.0,2.0,2.3,2.3)

  @test isa(mtv,MTensorValue{2,4})
  @test mtv[2,1] == 2.0

  p = Point(1.0,2.0,2.3)

  @test isa(p,Point{3})

  mp = MPoint(1.0,2.0,2.3)

  @test isa(mp,MPoint{3})

end

@testset "FieldValuesOuter" begin


end





end # module FieldValuesTests
