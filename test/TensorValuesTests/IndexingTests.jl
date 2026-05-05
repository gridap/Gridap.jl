module IndexingTests

using Test
using Gridap.TensorValues
using StaticArrays

@test Base.IndexStyle(MultiValue) == IndexCartesian()
@test Base.IndexStyle(VectorValue()) == IndexCartesian()

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
@test_throws BoundsError v[end+1]
@test_throws BoundsError v[0]

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
@test_throws BoundsError t[end+1]
@test_throws BoundsError t[0,0]

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
K = SkewSymTensorValue{2}(21)
t = TensorValue(convert(SMatrix{2,2,Int},s))
p = TensorValue(convert(SMatrix{2,2,Int},q))
r = TensorValue(convert(SMatrix{2,2,Int},K))

@test size(s) == size(q) == size(K) == (2,2)
@test length(s) == length(q) == length(K) == 4
@test lastindex(s) == lastindex(q) == lastindex(K) == length(s)
@test s[end] == 22
@test q[end] == -11
@test K[begin] == K[end] == 0
@test_throws BoundsError s[end+1]
@test_throws BoundsError s[0,0]
@test_throws BoundsError q[end+1]
@test_throws BoundsError q[0,0]
@test_throws BoundsError K[end+1]
@test_throws BoundsError K[0,0]

for (k,i) in enumerate(eachindex(t))
    @test s[i] == t[k]
end
for (k,i) in enumerate(eachindex(p))
    @test q[i] == p[k]
end
for (k,i) in enumerate(eachindex(r))
    @test K[i] == r[k]
end

@test s[2,1] == q[2,1] == 21
@test s[2] == q[2] == -K[2] == 21
@test s[3] == q[3] == K[3] == 21
@test q[1] == -q[4]

for (k,si) in enumerate(t)
  @test si == s[k]
end
for (k,qi) in enumerate(p)
  @test qi == q[k]
end
for (k,Ki) in enumerate(r)
  @test Ki == K[k]
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

# Voigt and Mandel notation
# SymTensorValue storage order: (1,1),(1,2),...,(1,D),(2,2),...,(D,D) (row-major upper triangle)
# Voigt order: diagonals first (1,1),(2,2),...,(D,D), then off-diagonal (1,2),(1,3),...

# Second-order SymTensorValue (D=2, L=3)
# a[1,1]=1, a[1,2]=a[2,1]=2, a[2,2]=3
a_sym2 = SymTensorValue(1.0, 2.0, 3.0)
v_a2 = to_voigt(a_sym2)
@test v_a2[1] == a_sym2[1,1]  # first diagonal
@test v_a2[2] == a_sym2[2,2]  # second diagonal
@test v_a2[3] == a_sym2[1,2]  # off-diagonal
@test from_voigt(v_a2) == a_sym2

m_a2 = to_mandel(a_sym2)
@test m_a2[1] ≈ a_sym2[1,1]              # diagonal unchanged
@test m_a2[2] ≈ a_sym2[2,2]              # diagonal unchanged
@test m_a2[3] ≈ sqrt(2) * a_sym2[1,2]   # off-diagonal scaled
@test from_mandel(m_a2) ≈ a_sym2

# Mandel isometry: to_mandel(a)⋅to_mandel(b) = inner(a,b)
b_sym2 = SymTensorValue(4.0, 5.0, 6.0)
@test to_mandel(a_sym2) ⋅ to_mandel(b_sym2) ≈ inner(a_sym2, b_sym2)

# Second-order SymTensorValue (D=3, L=6)
a_sym3 = SymTensorValue(1.0,2.0,3.0,4.0,5.0,6.0)
# storage: a[1,1]=1,a[1,2]=2,a[1,3]=3,a[2,2]=4,a[2,3]=5,a[3,3]=6
# Voigt order: [a11,a22,a33,a12,a13,a23] = [1,4,6,2,3,5]
@test to_voigt(a_sym3) == VectorValue(1.0,4.0,6.0,2.0,3.0,5.0)
@test from_voigt(to_voigt(a_sym3)) == a_sym3
@test from_mandel(to_mandel(a_sym3)) ≈ a_sym3

b_sym3 = SymTensorValue(7.0,8.0,9.0,10.0,11.0,12.0)
@test to_mandel(a_sym3) ⋅ to_mandel(b_sym3) ≈ inner(a_sym3, b_sym3)

# Type promotion: Mandel promotes T only when necessary (e.g. Int → Float64)
@test eltype(to_voigt(SymTensorValue(1, 2, 3))) == Int
@test eltype(to_mandel(SymTensorValue(1, 2, 3))) == Float64
@test eltype(to_mandel(SymTensorValue(1.0, 2.0, 3.0))) == Float64

# Fourth-order SymFourthOrderTensorValue (D=2, M=3)
C2 = one(SymFourthOrderTensorValue{2,Float64})
V2 = to_voigt(C2)
@test isa(V2, TensorValue{3,3,Float64})
# Voigt of identity: diagonal 1, off-diagonal 0 except V[3,3]=1/2
@test V2 ≈ TensorValue{3,3}(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,0.5)
@test from_voigt(V2) ≈ C2

M2 = to_mandel(C2)
@test isa(M2, TensorValue{3,3,Float64})
# Mandel of identity: the 4th-order identity maps to the 2nd-order identity matrix
@test M2 ≈ TensorValue{3,3}(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0)
@test from_mandel(M2) ≈ C2

# Fourth-order round-trip with non-trivial tensor (D=3)
C3 = one(SymFourthOrderTensorValue{3,Float64})
@test from_voigt(to_voigt(C3)) ≈ C3
@test from_mandel(to_mandel(C3)) ≈ C3

end # module IndexingTests
