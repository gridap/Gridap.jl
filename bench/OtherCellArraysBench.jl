
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

    struct DummyCellArray <: Numa.OtherCellArrayFromUnaryOp{Float64,2}
      a::Numa.OtherCellArray{Float64,1}
    end
    
    Numa.inputcellarray(self::DummyCellArray) = self.a
    
    Numa.computesize(self::DummyCellArray,asize) = (2,asize[1])
    
    function Numa.computevals!(self::DummyCellArray,a,v)
      @inbounds for i in 1:size(a,1)
        v[1,i] = a[i]
        v[2,i] = a[i]
      end
    end

  end)

  b = DummyCellArray(a)

  #(v,state) = iterate(b)
  #@code_warntype iterate(b)
  #v, astate = state
  #anext = iterate(Numa.inputcellarray(b),astate)
  #@code_warntype Numa.iteratekernel(b,anext,v)

  print("OtherCellArrayFromUnaryOp ->"); @time doloop(b)
  print("OtherCellArrayFromUnaryOp ->"); @time doloop(b)

end
