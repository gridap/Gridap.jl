
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

end
