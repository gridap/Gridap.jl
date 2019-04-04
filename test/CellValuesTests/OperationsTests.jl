
l = 10
sv = 1.0
scv = TestCellValue(sv,l)

sv2 = 1.1
scv2 = TestCellValue(sv2,l)

sa = [sv, sv, sv]
sca = TestCellValue(sa,l)

sa2 = [sv sv; sv sv; sv sv]
sca2 = TestCellValue(sa2,l)


@testset "Interfaces" begin

  @test eltype(scv) == Float64
  @test eltype(sca) == Array{Float64,1}
  @test eltype(sca2) == Array{Float64,2}

  s = string(sca)
  s0 = """
  1 -> [1.0, 1.0, 1.0]
  2 -> [1.0, 1.0, 1.0]
  3 -> [1.0, 1.0, 1.0]
  4 -> [1.0, 1.0, 1.0]
  5 -> [1.0, 1.0, 1.0]
  6 -> [1.0, 1.0, 1.0]
  7 -> [1.0, 1.0, 1.0]
  8 -> [1.0, 1.0, 1.0]
  9 -> [1.0, 1.0, 1.0]
  10 -> [1.0, 1.0, 1.0]
  """
  @test s == s0

  @test cellsize(sca2) == (3,2)
  @test cellsize(sca2,1) == 3
  @test cellsize(sca2,2) == 2
  @test celllength(sca2) == 6

end

using Numa.CellValues: CellValueFromUnaryOp
using Numa.CellValues: CellValueFromBinaryOp

@testset "Operations" begin

  for op in (:+,:-)
    @eval begin
      scv3 = $op(scv)
      @test isa(scv3,CellValueFromUnaryOp{Float64,typeof($op)})
      @test length(scv3) == l
      for vi in scv3
        @assert vi == $op(sv)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(inner),:(outer))
    @eval begin
      scv3 = $op(scv,scv2)
      @test isa(scv3,CellValueFromBinaryOp{Float64,typeof($op)})
      @test length(scv3) == l
      for vi in scv3
        @assert vi == $op(sv,sv2)
      end
    end
  end

end
