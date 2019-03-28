using LinearAlgebra

@time @testset "CellArrays" begin 

  using Numa.CellArrays
  using Numa.CellArrays: leftcellarray
  using Numa.CellArrays: rightcellarray

  aa = [1.0,2.0,2.1]
  aa2 = [0.0,2.1,1.1]
  bb = [aa';aa']
  l = 10
  a = ConstantCellArray(aa,l)
  a2 = ConstantCellArray(aa2,l)

  tv = TensorValue{2,4}(0.0,1.0,2.0,2.0)
  tt = [tv, tv, 4*tv, -1*tv]
  t = ConstantCellArray(tt,l)

  @testset "ConstantCellArray" begin

    @test length(a) == l
    @test cellsize(a) == size(aa)
    @test cellsize(a,1) == size(aa,1)
    @test celllength(a) == length(aa)
    for ar in a
      @assert ar == aa
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

  @testset "CachedArray" begin

    using Numa.CellArrays: CachedArray

    x = rand(10,5)

    y = CachedArray(x)

    @test y == x

    z = y + x

  end

  @testset "CellArrayFromUnaryOp" begin

    eval(quote


      struct DummyCellArray{C} <: CellArrayFromUnaryOp{C,Float64,2}
        a::C
      end

      import Numa.CellArrays: inputcellarray
      import Numa.CellArrays: computesize
      import Numa.CellArrays: computevals!
      
      inputcellarray(self::DummyCellArray) = self.a
      
      computesize(self::DummyCellArray,asize) = (2,asize[1])
      
      function computevals!(self::DummyCellArray,a,v)
        @inbounds for j in 1:size(a,1)
          for i in 1:2
            v[i,j] = a[j]
          end
        end
      end

    end)

    b = DummyCellArray(a)

    @test inputcellarray(b) === a
    @test length(b) == l
    @test cellsize(b) == (2,size(aa,1))
    @test cellsize(b,1) == 2
    @test cellsize(b,2) == size(aa,1)
    @test celllength(b) == 2*length(aa)
    for br in b
      @assert br == bb
    end

  end

  @testset "CellArrayFromDet" begin

    using Numa.CellArrays: CellArrayFromDet

    dett = [ det(tti) for tti in tt ]

    b = CellArrayFromDet{typeof(t),Float64,1}(t)

    @test inputcellarray(b) === t
    @test length(b) == l
    @test cellsize(b) == size(tt)
    @test cellsize(b,1) == size(tt,1)
    @test celllength(b) == size(tt,1)
    for br in b
      @assert br == dett
    end

    c = det(t)

    @test b == c

    @test isa(c,ConstantCellArray)

  end

  @testset "CellArrayFromInv" begin

    using Numa.CellArrays: CellArrayFromInv

    invt = [ inv(tti) for tti in tt ]

    b = CellArrayFromInv{typeof(t),typeof(tv),1}(t)

    @test inputcellarray(b) === t
    @test length(b) == l
    @test cellsize(b) == size(tt)
    @test cellsize(b,1) == size(tt,1)
    @test celllength(b) == size(tt,1)
    for br in b
      @assert br == invt
    end

    c = inv(t)

    @test b == c

    @test isa(c,ConstantCellArray)

  end

  @testset "CellArrayFromSum" begin

    using Numa.CellArrays: CellArrayFromSum

    b = CellArrayFromSum{typeof(a),typeof(a2),Float64,1}(a,a2)

    @test leftcellarray(b) === a
    @test rightcellarray(b) === a2
    @test length(b) == l
    @test cellsize(b) == size(aa)
    @test cellsize(b,1) == size(aa,1)
    @test celllength(b) == size(aa,1)
    for br in b
      @assert br == aa + aa2
    end

    c = a + a2

    @test b == c

    @test isa(c,ConstantCellArray)

  end

  @testset "CellArrayFromSum" begin

    using Numa.CellArrays: CellArrayFromSub

    b = CellArrayFromSub{typeof(a),typeof(a2),Float64,1}(a,a2)

    @test leftcellarray(b) === a
    @test rightcellarray(b) === a2
    @test length(b) == l
    @test cellsize(b) == size(aa)
    @test cellsize(b,1) == size(aa,1)
    @test celllength(b) == size(aa,1)
    for br in b
      @assert br == aa - aa2
    end

    c = a - a2

    @test b == c

    @test isa(c,ConstantCellArray)

  end

  @testset "CellArrayFromMul" begin

    using Numa.CellArrays: CellArrayFromMul

    b = CellArrayFromMul{typeof(a),typeof(a2),Float64,1}(a,a2)

    @test leftcellarray(b) === a
    @test rightcellarray(b) === a2
    @test length(b) == l
    @test cellsize(b) == size(aa)
    @test cellsize(b,1) == size(aa,1)
    @test celllength(b) == size(aa,1)
    for br in b
      @assert br == aa .* aa2
    end

    c = a * a2

    @test b == c

    @test isa(c,ConstantCellArray)

  end

  @testset "CellArrayFromDiv" begin

    using Numa.CellArrays: CellArrayFromDiv

    b = CellArrayFromDiv{typeof(a),typeof(a2),Float64,1}(a,a2)

    @test leftcellarray(b) === a
    @test rightcellarray(b) === a2
    @test length(b) == l
    @test cellsize(b) == size(aa)
    @test cellsize(b,1) == size(aa,1)
    @test celllength(b) == size(aa,1)
    for br in b
      @assert br == aa ./ aa2
    end

    c = a / a2

    @test b == c

    @test isa(c,ConstantCellArray)

  end

end
