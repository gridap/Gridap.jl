
let

  aa = [1.0,2.0,2.1]
  l = 10
  a = Numa.OtherConstantCellArray(aa,l)

  @test length(a) == l
  @test maxsize(a) == size(aa)

  for (ar,s) in a
    @assert s == size(ar)
    @assert ar == aa
  end

  cs(asize) = (2,asize[1])

  function cv!(avals,asize,vals,s)
    @inbounds for i in 1:asize[1]
      vals[1,i] = avals[i]
      vals[2,i] = avals[i]
    end
  end

  b = Numa.OtherCellArrayFromUnaryOp{Float64,1,Float64,2}(a,cv!,cs)

  bb = [aa';aa']

  @test length(b) == l
  @test maxsize(b) == size(bb)

  for (br,s) in b
    @assert s == size(br)
    @assert br == bb
  end

  c = Numa.DummyCellArray(a)

  @test length(c) == l
  @test maxsize(c) == size(bb)

  for (br,s) in c
    @assert s == size(br)
    @assert br == bb
  end

end
