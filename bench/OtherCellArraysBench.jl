
const N = 1000000

let

  println("+++ OtherCellArraysBench ( length = $N ) +++")

  function doloop(a)
    for ai in a
    end
  end

  aa = [1.0,2.0,2.1]
  a = Numa.OtherConstantCellArray(aa,N)

  print("OtherConstantCellArray ->"); @time doloop(a)
  print("OtherConstantCellArray ->"); @time doloop(a)

  eval(quote

    struct DummyCellArray{C} <: Numa.OtherCellArrayFromUnaryOp{C,Float64,2}
      a::C
    end
    
    Numa.inputcellarray(self::DummyCellArray) = self.a
    
    Numa.computesize(self::DummyCellArray,asize) = (2,asize[1])
    
    function Numa.computevals!(self::DummyCellArray,a,v)
      @inbounds for j in 1:size(a,1)
        for i in 1:2
          v[i,j] = a[j]
        end
      end
    end

  end)

  b = DummyCellArray(a)

  print("OtherCellArrayFromUnaryOp ->"); @time doloop(b)
  print("OtherCellArrayFromUnaryOp ->"); @time doloop(b)

end
