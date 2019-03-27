using LinearAlgebra

@time @testset "OtherCellArrays" begin 

  using Numa.OtherCellArrays

  aa = [1.0,2.0,2.1]
  bb = [aa';aa']
  l = 10
  a = OtherConstantCellArray(aa,l)

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


      struct DummyCellArray{C} <: OtherCellArrayFromUnaryOp{C,Float64,2}
        a::C
      end

      import Numa.OtherCellArrays: inputcellarray
      import Numa.OtherCellArrays: computesize
      import Numa.OtherCellArrays: computevals!
      
      inputcellarray(self::DummyCellArray) = self.a
      
      computesize(self::DummyCellArray,asize) = (2,asize[1])
      
      function computevals!(self::DummyCellArray,a,asize,v,vsize)
        @assert vsize == (2,asize[1])
        @inbounds for j in 1:asize[1]
          for i in 1:2
            v[i,j] = a[j]
          end
        end
      end

    end)

    b = DummyCellArray(a)

    @test inputcellarray(b) === a
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

  @testset "OtherCellArrayFromDet" begin

    using Numa.OtherCellArrays: OtherCellArrayFromDet

    tv = TensorValue{2,4}(0.0,1.0,2.0,2.0)
    tt = [tv, tv, 4*tv, -1*tv]
    dett = [ det(tti) for tti in tt ]
    t = OtherConstantCellArray(tt,l)

    b = OtherCellArrayFromDet{typeof(t),Float64,1}(t)

    @test inputcellarray(b) === t
    @test length(b) == l
    @test maxsize(b) == size(tt)
    @test maxsize(b,1) == size(tt,1)
    @test eltype(b) == Array{Float64,1}
    @test maxlength(b) == size(tt,1)
    for (br,brs) in b
      @assert br == dett
      @assert brs == size(tt)
    end

  end

  @testset "DetOfConstantCellArray" begin

    tv = TensorValue{2,4}(0.0,1.0,2.0,2.0)
    tt = [tv, tv, 4*tv, -1*tv]
    dett = [ det(tti) for tti in tt ]
    t = OtherConstantCellArray(tt,l)

    b = det(t)

    @test isa(b,OtherConstantCellArray)

    @test length(b) == l
    @test maxsize(b) == size(tt)
    @test maxsize(b,1) == size(tt,1)
    @test eltype(b) == Array{Float64,1}
    @test maxlength(b) == size(tt,1)
    for (br,brs) in b
      @assert br == dett
      @assert brs == size(tt)
    end

  end

end
