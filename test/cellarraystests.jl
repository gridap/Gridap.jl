
struct DummyCellArray{T,N} <: IndexableCellArray{T,N}
  a::Array{T,N}
  l::Int
end

Base.getindex(self::DummyCellArray,cell::Int) = self.a

Base.length(self::DummyCellArray) = self.l

a = [3.0,0.1,2.0]
l = 10

dca = DummyCellArray(a,l)

@test isa(dca,CellArray)

@test eltype(DummyCellArray{Float64,1}) == Array{Float64,1}

@test eltype(dca) == Array{Float64,1}

@test length(dca) == l

for b in dca
  @test b == a
end

#println(dca)

