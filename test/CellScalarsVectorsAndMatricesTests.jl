
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

