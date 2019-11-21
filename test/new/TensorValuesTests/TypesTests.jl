module TypesTests

using Gridap.TensorValues
using Test
using StaticArrays

# Constructors (MultiValue)

a = MArray{Tuple{3,2}}((1,2,3,4,5,6))
v = MultiValue(a)
@test isa(v,MultiValue{Tuple{3,2},Int})
@test v.array.data === a.data

a = SArray{Tuple{3,2}}((1,2,3,4,5,6))

v = MultiValue(a)
@test isa(v,MultiValue{Tuple{3,2},Int})
@test v.array === a

v = MultiValue{Tuple{3,2}}((1,2,3,4,5,6))
@test isa(v,MultiValue{Tuple{3,2},Int})
@test v.array == a

v = MultiValue{Tuple{3,2}}(1,2,3,4,5,6)
@test isa(v,MultiValue{Tuple{3,2},Int})
@test v.array == a

v = MultiValue{Tuple{3,2},Float64}((1,2,3,4,5,6))
@test isa(v,MultiValue{Tuple{3,2},Float64})
@test v.array == a

v = MultiValue{Tuple{3,2},Float64}(1,2,3,4,5,6)
@test isa(v,MultiValue{Tuple{3,2},Float64})
@test v.array == a

a = SVector(1)
v = MultiValue{Tuple{1}}((1,))
@test isa(v,MultiValue{Tuple{1},Int})
@test v.array == a

v = MultiValue{Tuple{1}}(1)
@test isa(v,MultiValue{Tuple{1},Int})
@test v.array == a

a = SMatrix{1,1}(1)
v = MultiValue{Tuple{1,1}}(1)
@test isa(v,MultiValue{Tuple{1,1},Int})
@test v.array == a

a = SVector{0,Int}()
v = MultiValue{Tuple{0},Int}(())
@test isa(v,MultiValue{Tuple{0},Int})
@test v.array == a

a = SMatrix{0,0,Int}()
v = MultiValue{Tuple{0,0},Int}()
@test isa(v,MultiValue{Tuple{0,0},Int})
@test v.array == a

# Constructors (TensorValue)

a = SMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,Int})
@test t.array == [1 3;2 4]

a = MMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,Int})
@test t.array == [1 3;2 4]

t = TensorValue{2}((1,2,3,4))
@test isa(t,TensorValue{2,Int})
@test t.array == [1 3;2 4]

t = TensorValue{2,Float64}((1,2,3,4))
@test isa(t,TensorValue{2,Float64})
@test t.array == [1 3;2 4]

t = TensorValue{2}(1,2,3,4)
@test isa(t,TensorValue{2,Int})
@test t.array == [1 3;2 4]

t = TensorValue{2,Float64}(1,2,3,4)
@test isa(t,TensorValue{2,Float64})
@test t.array == [1 3;2 4]

t = TensorValue(1,2,3,4)
@test isa(t,TensorValue{2,Int})
@test t.array == [1 3;2 4]

t = TensorValue((1,2,3,4))
@test isa(t,TensorValue{2,Int})
@test t.array == [1 3;2 4]

t = TensorValue{0,Int}()
@test isa(t,TensorValue{0,Int})
@test t.array == zeros(0,0)

t = TensorValue{1}(10)
@test isa(t,TensorValue{1,Int})
@test t.array == 10*ones(1,1)

t = TensorValue{1}((10,))
@test isa(t,TensorValue{1,Int})
@test t.array == 10*ones(1,1)

t = TensorValue{1,Float64}(10)
@test isa(t,TensorValue{1,Float64})
@test t.array == 10*ones(1,1)

t = TensorValue{1,Float64}((10,))
@test isa(t,TensorValue{1,Float64})
@test t.array == 10*ones(1,1)

# Constructors (VectorValue)

a = SVector(1)
g = VectorValue(a)
@test isa(g,VectorValue{1,Int})
@test g.array == [1,]

a = SVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

a = MVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

g = VectorValue{4}((1,2,3,4))
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

g = VectorValue{1}((1,))
@test isa(g,VectorValue{1,Int})
@test g.array == [1,]

g = VectorValue{0,Int}(())
@test isa(g,VectorValue{0,Int})
@test g.array == []

g = VectorValue{4}(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

g = VectorValue{1}(1)
@test isa(g,VectorValue{1,Int})
@test g.array == [1,]

g = VectorValue{1,Float64}(1)
@test isa(g,VectorValue{1,Float64})
@test g.array == [1,]

g = VectorValue{1,Float64}((1,))
@test isa(g,VectorValue{1,Float64})
@test g.array == [1,]

g = VectorValue{0,Int}()
@test isa(g,VectorValue{0,Int})
@test g.array == []

g = VectorValue{4,Float64}((1,2,3,4))
@test isa(g,VectorValue{4,Float64})
@test g.array == [1,2,3,4]

g = VectorValue{4}(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

g = VectorValue{4,Float64}(1,2,3,4)
@test isa(g,VectorValue{4,Float64})
@test g.array == [1,2,3,4]

g = VectorValue(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

g = VectorValue((1,2,3,4))
@test isa(g,VectorValue{4,Int})
@test g.array == [1,2,3,4]

g = VectorValue(1)
@test isa(g,VectorValue{1,Int})
@test g.array == [1,]

# Initializers

z = zero(MultiValue{Tuple{3,2},Int,2,6})
@test isa(z,MultiValue{Tuple{3,2},Int,2,6})
@test z.array == zeros(Int,(3,2))
s = zero(z)
@test s.array == zeros(Int,(3,2))

z = zero(TensorValue{3,Int,9})
@test isa(z,TensorValue{3,Int,9})
@test z.array == zeros(Int,(3,3))

z = zero(VectorValue{3,Int})
@test isa(z,VectorValue{3,Int})
@test z.array == zeros(Int,3)

z = one(TensorValue{3,Int,9})
@test isa(z,TensorValue{3,Int,9})
@test z.array == [1 0 0; 0 1 0; 0 0 1]
s = one(z)
@test s.array == [1 0 0; 0 1 0; 0 0 1]

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
V = MultiValue{Tuple{3,2},Int,2,6}
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

v = MultiValue{Tuple{3,2},Float64}(1,2,3,4,5,6)
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
@test change_eltype(a,Float64) == TensorValue{2,Float64,4}

@test change_eltype(1,Float64) == Float64

end # module TypesTests
