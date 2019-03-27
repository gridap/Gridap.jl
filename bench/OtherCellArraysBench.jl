
const N = 1000000

let

  using Numa.OtherCellArrays

  println("+++ OtherCellArraysBench ( length = $N ) +++")

  function doloop(a)
    for ai in a
    end
  end

  aa = [1.0,2.0,2.1]
  a = OtherConstantCellArray(aa,N)
  aa2 = [0.0,2.1,1.1]
  a2 = OtherConstantCellArray(aa2,N)

  print("OtherConstantCellArray ->"); @time doloop(a)
  print("OtherConstantCellArray ->"); @time doloop(a)

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
  t = OtherConstantCellArray(tt,N)
  c = Numa.OtherCellArrays.OtherCellArrayFromDet{typeof(t),Float64,1}(t)

  print("OtherCellArrayFromDet ->"); @time doloop(c)
  print("OtherCellArrayFromDet ->"); @time doloop(c)

  d = Numa.OtherCellArrays.OtherCellArrayFromInv{typeof(t),typeof(tv),1}(t)

  print("OtherCellArrayFromInv ->"); @time doloop(d)
  print("OtherCellArrayFromInv ->"); @time doloop(d)

  e = Numa.OtherCellArrays.OtherCellArrayFromSum{typeof(a),typeof(a2),Float64,1}(a,a2)

  print("OtherCellArrayFromSum ->"); @time doloop(e)
  print("OtherCellArrayFromSum ->"); @time doloop(e)

end
