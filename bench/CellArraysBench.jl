
const N = 1000000

let

  using Numa.CellArrays

  println("+++ CellArraysBench ( length = $N ) +++")

  function doloop(a)
    for ai in a
    end
  end

  aa = [1.0,2.0,2.1]
  a = ConstantCellArray(aa,N)
  aa2 = [0.0,2.1,1.1]
  a2 = ConstantCellArray(aa2,N)
  bb = [aa';aa']
  z = ConstantCellArray(bb,N)

  print("ConstantCellArray ->"); @time doloop(a)
  print("ConstantCellArray ->"); @time doloop(a)

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

  print("CellArrayFromUnaryOp ->"); @time doloop(b)
  print("CellArrayFromUnaryOp ->"); @time doloop(b)


  tv = TensorValue{2,4}(0.0,1.0,2.0,2.0)
  tt = [tv, tv, 4*tv, -1*tv]
  t = ConstantCellArray(tt,N)
  c = Numa.CellArrays.CellArrayFromDet{typeof(t),Float64,1}(t)

  print("CellArrayFromDet ->"); @time doloop(c)
  print("CellArrayFromDet ->"); @time doloop(c)

  d = Numa.CellArrays.CellArrayFromInv{typeof(t),typeof(tv),1}(t)

  print("CellArrayFromInv ->"); @time doloop(d)
  print("CellArrayFromInv ->"); @time doloop(d)

  e = Numa.CellArrays.CellArrayFromSum{typeof(a),typeof(a2),Float64,1}(a,a2)

  print("CellArrayFromSum ->"); @time doloop(e)
  print("CellArrayFromSum ->"); @time doloop(e)

  f = Numa.CellArrays.CellArrayFromCellSum{2,1,typeof(z),Float64}(z)

  print("CellArrayFromCellSum ->"); @time doloop(f)
  print("CellArrayFromCellSum ->"); @time doloop(f)

end
