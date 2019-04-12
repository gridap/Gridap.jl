
@testset "Wrappers" begin

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

  #       1,2,3,4,5,6,7,8,9,0,1,2
  data = [2,3,1,3,6,7,3,2,5,6,3,4]
  ptrs = [1,4,4,7,13]

  ca = CellVectorFromDataAndPtrs(data,ptrs)

  @test length(ca) == length(ptrs)-1
  @test ca[1] == data[1:3]
  @test ca[2] == data[4:3]
  @test ca[3] == data[4:6]
  @test ca[4] == data[7:12]

  @test cellsize(ca) == (6,)

  ca = CellVectorFromDataAndStride(data,3)

  @test length(ca) == 4
  @test ca[1] == data[1:3]
  @test ca[2] == data[4:6]
  @test ca[3] == data[7:9]
  @test ca[4] == data[10:12]

  @test cellsize(ca) == (3,)

end
