
let

  l = 10
  a = [1.1,2.1,2.4]
  b = [3.1,2.3,2.4]
  test = ConstantCellArray{Float64,1}(a,l)
  trial = ConstantCellArray{Float64,1}(b,l)
  
  values = inner(test,trial)
  
  @test length(values) == l
  @test maxsize(values) == (3,)
  
  for vals in values
    for (i,v) in enumerate(vals)
      @test v == a[i]*b[i]
    end
  end

end

let

  l = 10
  a = [1.1 2.1 2.4; 2.3 1.2 3.1; 3.2 4.3 6.5; 2.3 1.2 3.1]
  b = [3.1,2.3,2.4]
  test = ConstantCellArray{Float64,2}(a,l)
  trial = ConstantCellArray{Float64,1}(b,l)
  
  values = inner(test,trial)
  
  @test length(values) == l
  @test maxsize(values) == (4,3)
  
  for vals in values
    for i in 1:size(a,2)
      for j in 1:size(a,1)
        @test vals[j,i] == a[j,i]*b[i]
      end
    end
  end

end

