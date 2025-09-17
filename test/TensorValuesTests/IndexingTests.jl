module IndexingTests

using Test
using Gridap.TensorValues
using StaticArrays

a = (3,4,5,1)

v = VectorValue{4}(a)

@test IndexStyle(v) == IndexCartesian()
@test IndexStyle(typeof(v)) == IndexCartesian()

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

@test IndexStyle(t) == IndexCartesian()
@test IndexStyle(typeof(t)) == IndexCartesian()

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

# Index and size APIs
v = SMatrix{2,3}(1:6...)
w = TensorValue{2,3}(1:6...)
@test SArray(w) === v
@test w === MultiValue(v)

@test CartesianIndices(w) === CartesianIndices(v)
@test LinearIndices(w)    === LinearIndices(v)

@test eachindex(IndexCartesian(), w) === eachindex(IndexCartesian(), v)
@test eachindex(IndexLinear(), w)    === eachindex(IndexLinear(), v)

@test keys(w) === keys(v)
@test keys(IndexCartesian(), w) === keys(IndexCartesian(), v)
@test keys(IndexLinear(), w)    === keys(IndexLinear(), v)

@test lastindex(w)   === lastindex(v)
@test lastindex(w,1) === lastindex(v,1)
@test lastindex(w,2) === lastindex(v,2)

@test w[] === w # Number convention

# "non scalar" indexing

# Statically inferable index sets / ranges
@test w[:] === MultiValue(v[:])
@test w[:,:] === MultiValue(v[:,:])

@test w[begin,:] isa VectorValue{3}
@test w[begin,:] === MultiValue(v[begin,:])
@test w[:,end] === MultiValue(v[:,end])
@test w[end,end] == v[end,end]

So0 = SOneTo(0)
So2 = SOneTo(2)
@test w[So0,:] isa TensorValue{0,3}
@test w[So2,:]   === MultiValue(v[So2,:])
@test w[So2,So2] === MultiValue(v[So2,So2])
@test w[So0,So2] === MultiValue(v[So0,So2])

@test w[So0,1] isa VectorValue{0}
@test w[end,So0] === MultiValue(v[end,So0])
@test w[So2,begin] === MultiValue(v[So2,begin])

v0 = SVector{0,Int}()
v1 = SVector(1, 2)
v2 = SVector(1, 3)
@test w[v0,:] isa TensorValue{0,3}
@test w[v1,:]  == MultiValue(v[v1,:])
@test w[:,v2]  == MultiValue(v[:,v2])
@test w[v1,v2] == MultiValue(v[v1,v2])
@test w[v0,v2] == MultiValue(v[v0,v2])

@test w[v0,1] isa VectorValue{0}
@test w[end,v2] === MultiValue(v[end,v2])
@test w[v1,begin] === MultiValue(v[v1,begin])

# Non-statically inferable index sets / ranges

@test w[ [] ] == v[ [] ] == Int[] # vector
@test w[ [], [] ] == v[  [], []  ] == Int[;;] # matrix

o0 = Base.OneTo(0)
o2 = Base.OneTo(2)

@test w[o0,:] isa Matrix
@test w[o2,:]   == v[o2,:]
@test w[o2,o2] == v[o2,o2]
@test w[o0,o2] == v[o0,o2]

@test w[o0,1] isa Vector
@test w[end,o0] == v[end,o0]
@test w[o2,begin] == v[o2,begin]

v0 = Int[]
v1 = Int[1, 2]
v2 = Int[1, 3]
@test w[v0,:] isa Matrix{Int}
@test w[v1,:]  == v[v1,:]
@test w[:,v2]  == v[:,v2]
@test w[v1,v2] == v[v1,v2]
@test w[v0,v2] == v[v0,v2]

@test w[v0,1] isa Vector{Int}
@test w[end,v2] == v[end,v2]
@test w[v1,begin] == v[v1,begin]

mask = isodd.(v)
@test w[mask] == v[mask]
mask = Matrix(mask)
@test w[mask] == v[mask]

# Some BoundsError tests

# scalar indexing
t = ThirdOrderTensorValue(1:8...)
@test_throws BoundsError t[1,2]
@test_throws BoundsError t[CartesianIndex(1,2)]

# static
@test_throws BoundsError t[:,2]
@test_throws BoundsError t[v1,2]
@test_throws BoundsError t[v1,v1]

# dynamic
@test_throws BoundsError t[o2,2]
@test_throws BoundsError t[v1,2]
@test_throws BoundsError t[v1,v1]

end # module IndexingTests
