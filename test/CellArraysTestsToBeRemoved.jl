
let

  a = [3.0,0.1,2.0]
  l = 10
  
  dca = ConstantCellArray(a,l)
  
  @test isa(dca,CellArray)
  
  @test eltype(ConstantCellArray{Float64,1}) == Array{Float64,1}
  
  @test eltype(dca) == Array{Float64,1}
  
  @test length(dca) == l
  
  for b in dca
    @test b == a
  end
  
  s=string(dca)
  
  s_ref = """
  1 -> [3.0, 0.1, 2.0]
  2 -> [3.0, 0.1, 2.0]
  3 -> [3.0, 0.1, 2.0]
  4 -> [3.0, 0.1, 2.0]
  5 -> [3.0, 0.1, 2.0]
  6 -> [3.0, 0.1, 2.0]
  7 -> [3.0, 0.1, 2.0]
  8 -> [3.0, 0.1, 2.0]
  9 -> [3.0, 0.1, 2.0]
  10 -> [3.0, 0.1, 2.0]
  """
  
  @test s == s_ref

end

let

  l = 10
  aa = [3.0,0.1,2.0]
  bb = [1.0,0.0,2.0]
  
  a = ConstantCellArray(aa,l)
  b = ConstantCellArray(bb,l)

  T=Float64

  function computevals!(a,b,c)
    c .= a .+ b
    c
  end

  computesize(a,b) = a

  c = Numa.CellArrayFromBinaryOp{T,1,T,1,T,1}(a,b,computevals!,computesize)

  @test length(c) == l
  @test maxsize(c) == (3,)

  for ci in c
    @assert ci == aa + bb
  end


end

let

  l = 10
  aa = [3.0,0.1,2.0]
  
  a = ConstantCellArray(aa,l)

  T=Float64

  function computevals!(a,c)
    c .= 2a
    c
  end

  c = Numa.CellArrayFromUnaryOp(a,computevals!)

  @test length(c) == l
  @test maxsize(c) == (3,)

  for ci in c
    @assert ci == 2aa
  end


end
