
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


