
@testset "Wappers" begin

  v = [ Point{2}(i,i) for i in 1:10 ]

  cv = CellValueFromArray(v)

  @test length(cv) == length(v)
  @test size(cv) == size(v)
  for (cvi,vi) in zip(cv,v)
    @assert cvi === vi
  end
  @test IndexStyle(typeof(v)) == IndexStyle(typeof(cv))

  a = [ [i+1,i+2,i+3] for i in 1:10 ]

  ca = CellArrayFromArrayOfArrays(a)

  @test length(ca) == length(a)
  @test size(ca) == size(a)
  for (cai,ai) in zip(ca,a)
    @assert cai === ai
  end
  @test IndexStyle(typeof(a)) == IndexStyle(typeof(ca))

end

