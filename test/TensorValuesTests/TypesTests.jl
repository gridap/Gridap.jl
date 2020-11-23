module TypesTests

using Gridap.TensorValues
using Test
using Gridap.Arrays: get_array
using StaticArrays

# Constructors (TensorValue)

a = SMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

a = MMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

t = TensorValue{2}((1,2,3,4))
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

t = TensorValue{2}(1,2,3,4)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

t = TensorValue(1,2,3,4)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

t = TensorValue((1,2,3,4))
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

t = TensorValue{1}(10)
@test isa(t,TensorValue{1,1,Int})
@test convert(SMatrix{1,1,Int},t) == 10*ones(1,1)

t = TensorValue{1}((10,))
@test isa(t,TensorValue{1,1,Int})
@test convert(SMatrix{1,1,Int},t) == 10*ones(1,1)

t = TensorValue(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Float64})
@test convert(SMatrix{2,2,Float64},t) == [1 3;2 4]

t = TensorValue{2}(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Float64})
@test convert(SMatrix{2,2,Float64},t) == [1 3;2 4]

t = TensorValue{2,2}(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Float64})
@test convert(SMatrix{2,2,Float64},t) == [1 3;2 4]

t = TensorValue{2,2,Int}(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

# Constructors (SymTensorValue)

s = SymTensorValue( (11,21,22) )
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},s) == [11 21;21 22]

s = SymTensorValue(11,21,22)
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Float64},s) == [11.0 21.0;21.0 22.0]

s = SymTensorValue{2}( (11,21,22) )
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},s) == [11 21;21 22]

s = SymTensorValue{2}(11,21,22)
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Float64},s) == [11.0 21.0;21.0 22.0]

s = SymTensorValue{2,Int}( (11,21,22) )
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},s) == [11 21;21 22]

s = SymTensorValue{2,Float64}(11,21,22)
@test isa(s,SymTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},s) == [11.0 21.0;21.0 22.0]

s = SymTensorValue{0,Int}( () )
@test isa(s,SymTensorValue{0,Int})
@test convert(SMatrix{0,0,Int},s) == Array{Any,2}(undef,0,0)

s = SymTensorValue{0,Int}()
@test isa(s,SymTensorValue{0,Int})
@test convert(SMatrix{0,0,Int},s) == Array{Any,2}(undef,0,0)

s = SymTensorValue(11,21.0,22)
@test isa(s,SymTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},s) == [11.0 21.0;21.0 22.0]

s = SymTensorValue{2}(11,21.0,22)
@test isa(s,SymTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},s) == [11.0 21.0;21.0 22.0]

s = SymTensorValue{2,Int}(11,21.0,22)
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},s) == [11.0 21.0;21.0 22.0]

# Constructors (SymFourthOrderTensorValue)

s = SymFourthOrderTensorValue( (1111,1121,1122, 2111,2121,2122, 2211,2221,2222) )
@test isa(s,SymFourthOrderTensorValue{2,Int})
@test Tuple(s) == (1111,1121,1122, 2111,2121,2122, 2211,2221,2222)

s = SymFourthOrderTensorValue(1111,2111,2211, 1121,2121,2221, 1122,2122,2222)
@test isa(s,SymFourthOrderTensorValue{2,Int})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

s = SymFourthOrderTensorValue{2}( (1111,2111,2211, 1121,2121,2221, 1122,2122,2222) )
@test isa(s,SymFourthOrderTensorValue{2,Int})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

s = SymFourthOrderTensorValue{2}(1111,2111,2211, 1121,2121,2221, 1122,2122,2222)
@test isa(s,SymFourthOrderTensorValue{2,Int})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

s = SymFourthOrderTensorValue{2,Int}( (1111,2111,2211, 1121,2121,2221, 1122,2122,2222) )
@test isa(s,SymFourthOrderTensorValue{2,Int})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

s = SymFourthOrderTensorValue{2,Float64}(1111,2111,2211, 1121,2121,2221, 1122,2122,2222)
@test isa(s,SymFourthOrderTensorValue{2,Float64})
@test Tuple(s) == (1111.0,2111.0,2211.0, 1121.0,2121.0,2221.0, 1122.0,2122.0,2222.0)

s = SymFourthOrderTensorValue{0,Int}( () )
@test isa(s,SymFourthOrderTensorValue{0,Int})
@test Tuple(s) == ()

s = SymFourthOrderTensorValue{0,Int}()
@test isa(s,SymFourthOrderTensorValue{0,Int})
@test Tuple(s) == ()

s = SymFourthOrderTensorValue(1111,2111,2211.0, 1121,2121.0,2221, 1122,2122,2222)
@test isa(s,SymFourthOrderTensorValue{2,Float64})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

s = SymFourthOrderTensorValue{2}(1111,2111,2211.0, 1121,2121.0,2221, 1122,2122,2222)
@test isa(s,SymFourthOrderTensorValue{2,Float64})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

s = SymFourthOrderTensorValue{2,Int}(1111,2111,2211.0, 1121,2121.0,2221, 1122,2122,2222)
@test isa(s,SymFourthOrderTensorValue{2,Int})
@test Tuple(s) == (1111,2111,2211, 1121,2121,2221, 1122,2122,2222  )

# Constructors (VectorValue)

a = SVector(1)
g = VectorValue(a)
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1,Int},g) == [1,]

a = SVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

a = MVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test convert(MVector{4,Int},g) == [1,2,3,4]

g = VectorValue{4}((1,2,3,4))
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

g = VectorValue{1}((1,))
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1,Int},g) == [1,]

g = VectorValue{0,Int}(())
@test isa(g,VectorValue{0,Int})
@test convert(SVector{0,Int},g) == []

g = VectorValue{4}(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

g = VectorValue{1}(1)
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1,Int},g) == [1,]

g = VectorValue{1,Float64}(1)
@test isa(g,VectorValue{1,Float64})
@test convert(SVector{1,Float64},g) == [1,]

g = VectorValue{0,Int}()
@test isa(g,VectorValue{0,Int})
@test convert(SVector{0,Int},g) == []

g = VectorValue{4,Float64}((1,2,3,4))
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]

g = VectorValue{4,Float64}(1,2,3,4)
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]

g = VectorValue(1,2,3,4)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

g = VectorValue((1,2,3,4))
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

g = VectorValue(1)
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1,Int},g) == [1,]

g = VectorValue(1.0,2,3.0,4)
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]

g = VectorValue{4}(1.0,2,3.0,4)
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]

g = VectorValue{4,Int}(1.0,2,3.0,4)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

g = VectorValue((1.0,2,3.0,4))
@test isa(g,VectorValue{4,Float64})
@test convert(SVector{4,Float64},g) == [1,2,3,4]


# Initializers

z = zero(TensorValue{3,3,Int,9})
@test isa(z,TensorValue{3,3,Int,9})
@test convert(SMatrix{3,3,Int},z) == zeros(Int,(3,3))

z = zero(SymTensorValue{3,Int})
@test isa(z,SymTensorValue{3,Int,6})
@test convert(SMatrix{3,3,Int},z) == zeros(Int,(3,3))

z = zero(ThirdOrderTensorValue{3,3,3,Int,27})
@test isa(z,ThirdOrderTensorValue{3,3,3,Int,27})
@test Tuple(z) == Tuple(zeros(Int,(27)))

z = zero(SymFourthOrderTensorValue{2,Int})
@test isa(z,SymFourthOrderTensorValue{2,Int,9})
@test Tuple(z) == Tuple(zeros(Int,(9)))

z = zero(VectorValue{3,Int})
@test isa(z,VectorValue{3,Int})
@test convert(SVector{3,Int},z) == zeros(Int,3)

z = one(TensorValue{3,3,Int,9})
@test isa(z,TensorValue{3,3,Int,9})
@test convert(SMatrix{3,3,Int},z) == [1 0 0; 0 1 0; 0 0 1]
s = one(z)
@test convert(SMatrix{3,3,Int},s) == [1 0 0; 0 1 0; 0 0 1]

z = one(SymFourthOrderTensorValue{2,Int})
@test isa(z,SymFourthOrderTensorValue{2})
@test Tuple(z) == (1.,0.,0., 0.,0.5,0., 0.,0.,1.)

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

a = (11,21,22)
V = SymTensorValue{2,Int,3}
b = convert(V,a)
@test isa(b,V)
b = V[a,a,a,]
@test isa(b,Vector{V})

a = (1111,1121,1122, 2111,2121,2122, 2211,2221,2222)
V = SymFourthOrderTensorValue{2,Int,9}
b = convert(V,a)
@test isa(b,V)
@test b[2,2,2,1] == 2221
@test b[2,2,1,2] == 2221
@test b[2,1,2,1] == 2121
@test b[1,2,2,1] == 2121
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

v = SymTensorValue{3,Int}(1, 0, 0, 1, 0, 1)
s = "(1, 0, 0, 1, 0, 1)"
@test string(v) == s

v = SymFourthOrderTensorValue{2,Int}(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)
s = "(1111, 1121, 1122, 2111, 2121, 2122, 2211, 2221, 2222)"
@test string(v) == s

# Third order tensors

a = SArray{Tuple{2,2,2}}(1,2,3,4,5,6,7,8)
t = ThirdOrderTensorValue(a)
@test isa(t,ThirdOrderTensorValue{2,2,2,Int})

t = ThirdOrderTensorValue(1,2,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Int})

t = ThirdOrderTensorValue{2}(1,2,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Int})

t = ThirdOrderTensorValue{2,2,2}(1,2,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Int})

t = ThirdOrderTensorValue{2,2,2,Float64}(1,2,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Float64})

t = ThirdOrderTensorValue(1,2.0,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Float64})

t = ThirdOrderTensorValue{2}(1,2.0,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Float64})

t = ThirdOrderTensorValue{2,2,2}(1,2.0,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Float64})

t = ThirdOrderTensorValue{2,2,2,Int}(1,2.0,3,4,5,6,7,8)
@test isa(t,ThirdOrderTensorValue{2,2,2,Int})

# Misc

v = VectorValue(3,2,1)
m = mutable(v)
@test m == get_array(v)
@test isa(m,MVector)

v = TensorValue{2,3}(1,2,3,4,5,6)
m = mutable(v)
@test m == get_array(v)
@test isa(m,MMatrix)

M = Mutable(VectorValue{3,Int})
@test M == MVector{3,Int}
m = zero(M)
v = VectorValue(m)
@test isa(v,VectorValue{3,Int})


@test num_components(Int) == 1
@test num_components(Float64) == 1
@test num_components(1.0) == 1
@test num_components(1) == 1
@test num_components(VectorValue{3,Float64}) == 3
@test num_components(VectorValue(1,2,3)) == 3
@test num_components(TensorValue(1,2,3,4)) == 4
@test num_components(SymTensorValue(1,2,3)) == 4
@test num_components(SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)) == 16

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

a = SymTensorValue(11,21,22)
@test change_eltype(a,Float64) == SymTensorValue{2,Float64,3}
@test isa(Tuple(a),Tuple)
@test Tuple(a) == a.data

a = SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)
@test change_eltype(a,Float64) == SymFourthOrderTensorValue{2,Float64,9}
@test isa(Tuple(a),Tuple)
@test Tuple(a) == a.data

end # module TypesTests
