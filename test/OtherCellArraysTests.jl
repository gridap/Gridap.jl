
@time @testset "OtherCellArrays" begin 

  aa = [1.0,2.0,2.1]
  l = 10
  a = Numa.OtherConstantCellArray(aa,l)

  @testset "OtherConstantCellArray" begin

    @test length(a) == l
    @test maxsize(a) == size(aa)
    @test maxsize(a,1) == size(aa,1)
    @test eltype(a) == Array{Float64,1}
    @test maxlength(a) == length(aa)
    for (ar,ars) in a
      @assert ar == aa
      @assert ars == size(aa)
    end
    s = string(a)
    s0 = """
    1 -> [1.0, 2.0, 2.1]
    2 -> [1.0, 2.0, 2.1]
    3 -> [1.0, 2.0, 2.1]
    4 -> [1.0, 2.0, 2.1]
    5 -> [1.0, 2.0, 2.1]
    6 -> [1.0, 2.0, 2.1]
    7 -> [1.0, 2.0, 2.1]
    8 -> [1.0, 2.0, 2.1]
    9 -> [1.0, 2.0, 2.1]
    10 -> [1.0, 2.0, 2.1]
    """
    @test s == s0

  end

  @testset "OtherCellArrayFromUnaryOp" begin

    eval(quote


      struct DummyCellArray{C} <: Numa.OtherCellArrayFromUnaryOp{C,Float64,2}
        a::C
      end
      
      Numa.inputcellarray(self::DummyCellArray) = self.a
      
      Numa.computesize(self::DummyCellArray,asize) = (2,asize[1])
      
      function Numa.computevals!(self::DummyCellArray,a,asize,v,vsize)
        @assert vsize == (2,asize[1])
        @inbounds for j in 1:asize[1]
          for i in 1:2
            v[i,j] = a[j]
          end
        end
      end

    end)

    bb = [aa';aa']
    b = DummyCellArray(a)

    @test Numa.inputcellarray(b) === a
    @test length(b) == l
    @test maxsize(b) == (2,size(aa,1))
    @test maxsize(b,1) == 2
    @test maxsize(b,2) == size(aa,1)
    @test eltype(b) == Array{Float64,2}
    @test maxlength(b) == 2*length(aa)
    for (br,brs) in b
      @assert br == bb
      @assert brs == size(bb)
    end

  end

  #cs(asize) = (2,asize[1])

  #function cv!(avals,asize,vals,s)
  #  @inbounds for i in 1:asize[1]
  #    vals[1,i] = avals[i]
  #    vals[2,i] = avals[i]
  #  end
  #end

  #b = Numa.OtherCellArrayFromUnaryOp{Float64,1,Float64,2}(a,cv!,cs)

  #bb = [aa';aa']

  #@test length(b) == l
  #@test maxsize(b) == size(bb)

  #for (br,s) in b
  #  @assert s == size(br)
  #  @assert br == bb
  #end

  #c = Numa.DummyCellArray(a)

  #@test length(c) == l
  #@test maxsize(c) == size(bb)

  #for (br,s) in c
  #  @assert s == size(br)
  #  @assert br == bb
  #end

end
