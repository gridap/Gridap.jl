
l = 10
sv = 1.0
scv = TestCellValue(sv,l)

sv2 = 1.1
scv2 = TestCellValue(sv2,l)

sa = [sv, sv, sv]
sca = TestCellArray(sa,l)

sa2 = [sv sv; sv sv; sv sv]
sca2 = TestCellArray(sa2,l)


@testset "Interfaces" begin

  @test eltype(scv) == Float64

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
using Numa.CellValues: CellArrayFromBroadcastUnaryOp
using Numa.CellValues: CellArrayFromCellSum
using Numa.CellValues: CellArrayFromCellNewAxis
using Numa.CellValues: CellValueFromCellArrayReduce
using Numa.CellValues: CellArrayFromBoradcastBinaryOp

@testset "Operations" begin

  for op in (:+,:-,:(inv),:(det))
    @eval begin
      scv3 = $op(scv)
      @test isa(scv3,CellValueFromUnaryOp{Float64,typeof($op)})
      @test length(scv3) == l
      for vi in scv3
        @assert vi == $op(sv)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      scv3 = $op(scv,scv2)
      @test isa(scv3,CellValueFromBinaryOp{Float64,typeof($op)})
      @test length(scv3) == l
      for vi in scv3
        @assert vi == $op(sv,sv2)
      end
    end
  end

  for op in (:+,:-,:(inv),:(det))
    @eval begin
      sca3 = $op(sca)
      @test isa(sca3,CellArrayFromBroadcastUnaryOp{typeof($op)})
      @test length(sca3) == l
      for vi in sca3
        @assert vi == $op.(sa)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      sca3 = $op(sca,sca2)
      @test isa(sca3,CellArrayFromBoradcastBinaryOp{typeof($op),Float64,2})
      @test length(sca3) == l
      @test cellsize(sca3) == (3,2)
      for vi in sca3
        @assert vi == $op.(sa,sa2)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      sca3 = $op(scv,sca2)
      @test isa(sca3,CellArrayFromBoradcastBinaryOp{typeof($op),Float64,2})
      @test length(sca3) == l
      @test cellsize(sca3) == (3,2)
      for vi in sca3
        @assert vi == $op.(sv,sa2)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      sca3 = $op(sca2,scv)
      @test isa(sca3,CellArrayFromBoradcastBinaryOp{typeof($op),Float64,2})
      @test length(sca3) == l
      @test cellsize(sca3) == (3,2)
      for vi in sca3
        @assert vi == $op.(sa2,sv)
      end
    end
  end

  sca3 = cellsum(sca2,dim=2)
  @test isa(sca3,CellArrayFromCellSum{2,1,typeof(sca2),Float64})
  @test length(sca3) == l
  @test cellsize(sca3) == (3,)
  for vi in sca3
    @assert vi == reshape(sum(sa2,dims=2),(3,))
  end

  sca3 = cellsum(sca,dim=1)
  @test isa(sca3,CellValueFromCellArrayReduce{Float64,typeof(sum),typeof(sca)})
  @test length(sca3) == l
  for vi in sca3
    @assert vi == sum(sa)
  end

  sca3 = cellnewaxis(sca2,dim=2)
  @test isa(sca3,CellArrayFromCellNewAxis{2,typeof(sca2),Float64,3})
  @test length(sca3) == l
  @test cellsize(sca3) == (3,1,2)
  for vi in sca3
    @assert vi == reshape(sa2,(3,1,2))
  end

end
