module IndexingTests

using Test
using Gridap.TensorValues
using StaticArrays

a = (3,4,5,1)

v = VectorValue{4}(a)

@test eltype(v) == Int
@test eltype(typeof(v)) == Int

@test size(v) == (4,)
@test length(v) == 4

for (k,i) in enumerate(eachindex(v))
  @test v[i] == a[k]
end

t = TensorValue{2}(a)

@test size(t) == (2,2)
@test length(t) == 4

for (k,i) in enumerate(eachindex(t))
  @test t[i] == a[k]
end

@test t[2,1] == 4

@test t[2] == 4

for (k,ti) in enumerate(t)
  @test ti == a[k]
end

v = @SMatrix zeros(2,3)
w = MultiValue(v)
@test CartesianIndices(w) == CartesianIndices(v)
@test LinearIndices(w) == LinearIndices(v)

end # module IndexingTests
