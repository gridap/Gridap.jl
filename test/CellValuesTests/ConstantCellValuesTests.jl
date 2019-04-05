
l = 10
sv = 1.0
scv = ConstantCellValue(sv,l)

sv2 = 1.1
scv2 = ConstantCellValue(sv2,l)

sa = [sv, sv, sv]
sca = ConstantCellArray(sa,l)

sa2 = [sv sv; sv sv; sv sv]
sca2 = ConstantCellArray(sa2,l)

@testset "ConstantCellValues" begin

  @test isa(scv,CellValue{Float64})
  @test isa(scv,CellValue{Float64})
  @test length(scv) == l
  for vi in scv
    @assert vi === sv
  end
  @test scv == scv

  @test isa(sca,CellVector{Float64})
  @test isa(sca,CellArray{Float64,1})
  @test length(sca) == l
  @test cellsize(sca) == (3,)
  for vi in sca
    @assert vi === sa
  end
  @test sca == sca

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      scv3 = $op(scv,scv2)
      @test isa(scv3,ConstantCellValue{Float64})
      @test length(scv3) == l
      for vi in scv3
        @assert vi == $op(sv,sv2)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      sca3 = $op(sca,sca2)
      @test isa(sca3,ConstantCellArray{Float64,2})
      @test length(sca3) == l
      @test cellsize(sca3) == (3,2)
      for vi in sca3
        @assert vi == broadcast($op,sa,sa2)
      end
    end
  end

  for op in (:+,:-,:*,:/,:(outer))
    @eval begin
      sca3 = $op(scv,sca2)
      @test isa(sca3,ConstantCellArray{Float64,2})
      @test length(sca3) == l
      @test cellsize(sca3) == (3,2)
      for vi in sca3
        @assert vi == broadcast($op,sv,sa2)
      end
    end
  end

  for op in (:+,:-,:(det),:(inv))
    @eval begin
      sca3 = $op(sca2)
      @test isa(sca3,ConstantCellArray{Float64,2})
      @test length(sca3) == l
      @test cellsize(sca3) == (3,2)
      for vi in sca3
        @assert vi == broadcast($op,sa2)
      end
    end
  end

  sca3 = cellsum(sca2,dim=2)
  @test isa(sca3,ConstantCellArray{Float64,1})
  @test length(sca3) == l
  @test cellsize(sca3) == (3,)
  for vi in sca3
    @assert vi == reshape(sum(sa2,dims=2),(3,))
  end

  sca3 = cellsum(sca,dim=1)
  @test isa(sca3,ConstantCellValue{Float64})
  @test length(sca3) == l
  for vi in sca3
    @assert vi == sum(sa)
  end

  sca3 = cellnewaxis(sca2,dim=2)
  @test isa(sca3,ConstantCellArray{Float64,3})
  @test length(sca3) == l
  @test cellsize(sca3) == (3,1,2)
  for vi in sca3
    @assert vi == reshape(sa2,(3,1,2))
  end

end
