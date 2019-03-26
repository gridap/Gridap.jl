
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
    
    function Numa.computevals!(self::DummyCellArray,a,asize,v,vsize)
      @inbounds for j in 1:asize[1]
        for i in 1:2
          v[i,j] = a[j]
        end
      end
    end

  end)

  b = DummyCellArray(a)

  print("OtherCellArrayFromUnaryOp ->"); @time doloop(b)
  print("OtherCellArrayFromUnaryOp ->"); @time doloop(b)


  tv = TensorValue{2,4}(0.0,1.0,2.0,2.0)
  tt = [tv, tv, 4*tv, -1*tv]
  t = Numa.OtherConstantCellArray(tt,N)
  c = Numa.OtherConstantCellArrayFromDet{typeof(t),Float64,1}(t)

  print("OtherConstantCellArrayFromDet ->"); @time doloop(c)
  print("OtherConstantCellArrayFromDet ->"); @time doloop(c)


end
