module TypesTests

using Gridap.TensorValues
using Test
using StaticArrays

# Constructors (TensorValue)

a = SMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2},t) == [1 3;2 4]

a = MMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix,t) == [1 3;2 4]

t = TensorValue{2}((1,2,3,4))
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2},t) == [1 3;2 4]

t = TensorValue{2}(1,2,3,4)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2},t) == [1 3;2 4]

t = TensorValue(1,2,3,4)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2},t) == [1 3;2 4]

t = TensorValue((1,2,3,4))
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2},t) == [1 3;2 4]

t = TensorValue{1}(10)
@test isa(t,TensorValue{1,1,Int})
@test convert(SMatrix{1,1},t) == 10*ones(1,1)

t = TensorValue{1}((10,))
@test isa(t,TensorValue{1,1,Int})
@test convert(SMatrix,t) == 10*ones(1,1)

# Constructors (VectorValue)

a = SVector(1)
g = VectorValue(a)
@test isa(g,VectorValue{1,Int})
@test convert(SVector,g) == [1,]

a = SVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test convert(SVector,g) == [1,2,3,4]

a = MVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test convert(MVector,g) == [1,2,3,4]

g = VectorValue{4}((1,2,3,4))
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4},g) == [1,2,3,4]

g = VectorValue{1}((1,))
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1},g) == [1,]

g = VectorValue{0,Int}(())
@test isa(g,VectorValue{0,Int})
@test convert(SVector{0,Int},g) == []

g = VectorValue{4}(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4},g) == [1,2,3,4]

g = VectorValue{1}(1)
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1},g) == [1,]

g = VectorValue{1,Float64}(1)
@test isa(g,VectorValue{1,Float64})
@test convert(SVector{1,Float64},g) == [1,]

g = VectorValue{0,Int}()
@test isa(g,VectorValue{0,Int})
@test convert(SVector{0,Int},g) == []

g = VectorValue{4,Float64}((1,2,3,4))
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]

g = VectorValue{4}(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4},g) == [1,2,3,4]

g = VectorValue{4,Float64}(1,2,3,4)
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]

g = VectorValue(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test convert(SVector,g) == [1,2,3,4]

g = VectorValue((1,2,3,4))
@test isa(g,VectorValue{4,Int})
@test convert(SVector,g) == [1,2,3,4]

g = VectorValue(1)
@test isa(g,VectorValue{1,Int})
@test convert(SVector,g) == [1,]

# Initializers

z = zero(TensorValue{3,3,Int,9})
@test isa(z,TensorValue{3,3,Int,9})
@test convert(SMatrix,z) == zeros(Int,(3,3))

z = zero(VectorValue{3,Int})
@test isa(z,VectorValue{3,Int})
@test convert(SVector,z) == zeros(Int,3)

z = one(TensorValue{3,3,Int,9})
@test isa(z,TensorValue{3,3,Int,9})
@test convert(SMatrix,z) == [1 0 0; 0 1 0; 0 0 1]
s = one(z)
@test convert(SMatrix,s) == [1 0 0; 0 1 0; 0 0 1]

# Conversions

a = @SVector ones(Int,3)
b = convert(VectorValue{3,Int},a)
@test isa(b,VectorValue{3,Int})

a = ones(Int,3)
b = convert(VectorValue{3,Int},a)
@test isa(b,VectorValue{3,Int})

a = ones(Int,1)
b = convert(VectorValue{1,Int},a)
@test isa(b,VectorValue{1,Int})

a = (1,2,2,1,3,2)
V = TensorValue{3,2,Int,6}
b = convert(V,a)
@test isa(b,V)
b = V[a,a,a,]
@test isa(b,Vector{V})

# Misc operations on the type itself

V = VectorValue{3,Int}
@test length(V) == 3
@test size(V) == (3,)

V = VectorValue{3}
@test length(V) == 3
@test size(V) == (3,)

# Custom type printing

v = TensorValue{3,2,Float64}(1,2,3,4,5,6)
s = "(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)"
@test string(v) == s

# Misc

M = mutable(VectorValue{3,Int})
@test M == MVector{3,Int}
m = zero(M)
v = VectorValue(m)
@test isa(v,VectorValue{3,Int})

@test n_components(Int) == 1
@test n_components(Float64) == 1
@test n_components(1.0) == 1
@test n_components(1) == 1
@test n_components(VectorValue{3,Float64}) == 3
@test n_components(VectorValue(1,2,3)) == 3

a = VectorValue(1,2,3,4)
@test change_eltype(a,Float64) == VectorValue{4,Float64}

a = TensorValue(1,2,3,4)
@test change_eltype(a,Float64) == TensorValue{2,2,Float64,4}

@test change_eltype(1,Float64) == Float64

a = TensorValue(1,2,3,4)
@test isa(Tuple(a),Tuple)
@test Tuple(a) == a.data

p = VectorValue(2,3)
t = diagonal_tensor(p)
@test t == TensorValue(2,0,0,3)

p = VectorValue(1,2,3)
t = diagonal_tensor(p)
@test t == TensorValue(1,0,0,0,2,0,0,0,3)

end # module TypesTests
