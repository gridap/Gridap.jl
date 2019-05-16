module FieldValuesTests

using Test
using Gridap
using Gridap.FieldValues



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

  @test outer(1.2,3.4) == 1.2*3.4
  @test outer(Float64,Float64) == Float64

  vv = VectorValue(1.0,2.0,2.3)
  @test outer(1.3,vv) == 1.3*vv
  @test outer(Float64,typeof(vv)) == typeof(vv)

  tv = TensorValue(1.0,2.0,2.3,2.3)
  @test outer(1.3,tv) == 1.3*tv
  @test outer(Float64,typeof(tv)) == typeof(tv)

  vva = VectorValue(1.3,2.0,5.3)
  vvb = VectorValue(1.0,4.0,2.3)
  tvc = outer(vva,vvb)
  @test isa(tvc,TensorValue{3,9})
  @test outer(typeof(vva),typeof(vvb)) == typeof(tvc)
  for j in 1:3
    for i in 1:3
      @assert vva[i]*vvb[j] == tvc[i,j]
    end
  end

end

@testset "FieldValuesInner" begin

  @test inner(1.2,3.4) == 1.2*3.4
  @test inner(Float64,Float64) == Float64

  vva = VectorValue(1.3,2.0,5.3)
  vvb = VectorValue(1.0,4.0,2.3)
  c = inner(vva,vvb)
  @test isa(c,Float64)
  @test inner(typeof(vva),typeof(vvb)) == Float64
  @test c == vva[1]*vvb[1] + vva[2]*vvb[2] + vva[3]*vvb[3]

  tva = TensorValue(1.0,2.0,2.3,2.3)
  tvb = TensorValue(1.1,3.4,5.3,9.3)
  d = inner(tva,tvb)
  @test isa(d,Float64)
  @test inner(typeof(tva),typeof(tvb)) == Float64
  @test d == sum([ tva[i,j]*tvb[i,j] for j in 1:2 for i in 1:2 ])

end

end # module FieldValuesTests
