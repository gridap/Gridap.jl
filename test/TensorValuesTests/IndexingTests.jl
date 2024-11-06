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
@test lastindex(v) == length(v)
@test v[end] == a[end]

for (k,i) in enumerate(eachindex(v))
  @test v[i] == a[k]
end

t = TensorValue{2}(a)

@test size(t) == (2,2)
@test length(t) == 4
@test lastindex(t) == length(t)
@test t[end] == a[end]

for (k,i) in enumerate(eachindex(t))
  @test t[i] == a[k]
end

@test t[2,1] == 4

@test t[2] == 4

for (k,ti) in enumerate(t)
  @test ti == a[k]
end

s = SymTensorValue{2}(11,21,22)
q = SymTracelessTensorValue{2}(11,21)
t = TensorValue(convert(SMatrix{2,2,Int},s))
p = TensorValue(convert(SMatrix{2,2,Int},q))

@test size(s) == size(q) == (2,2)
@test length(s) == length(q) == 4
@test lastindex(s) == lastindex(q) == length(s)
@test s[end] == 22
@test q[end] == -11

for (k,i) in enumerate(eachindex(t))
    @test s[i] == t[k]
end
for (k,i) in enumerate(eachindex(p))
    @test q[i] == p[k]
end

@test s[2,1] == q[2,1] == 21
@test s[2] == q[2] == 21
@test q[1] == -q[4]

for (k,si) in enumerate(t)
  @test si == s[k]
end
for (k,qi) in enumerate(p)
  @test qi == q[k]
end

v = @SMatrix zeros(2,3)
w = TensorValue(v)
@test CartesianIndices(w) == CartesianIndices(v)
@test LinearIndices(w) == LinearIndices(v)

end # module IndexingTests
