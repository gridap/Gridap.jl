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
@test_throws DimensionMismatch convert(SMatrix{1,2,Int},t)
@test convert(SMatrix{2,2,Int},t) === SArray(t)
@test convert(MMatrix{2,2,Int},t) == MArray(t)
@test t === MultiValue(a)

a = MMatrix{2,2}(1,2,3,4)
t = TensorValue(a)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]
@test t === MultiValue(a)

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
@test convert(SMatrix{2,2,Float64},t) ===SMatrix(t)
@test convert(MMatrix{2,2,Float64},t) == MMatrix(t)

t = TensorValue{2}(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Float64})
@test convert(SMatrix{2,2,Float64},t) == [1 3;2 4]

t = TensorValue{2,2}(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Float64})
@test convert(SMatrix{2,2,Float64},t) == [1 3;2 4]

t = TensorValue{2,2,Int}(1,2.0,3,4)
@test isa(t,TensorValue{2,2,Int})
@test convert(SMatrix{2,2,Int},t) == [1 3;2 4]

a = [11.0 21.0; 12.0 22.0]
@test isa(a,AbstractArray{Float64,2})
t = convert(TensorValue{2,2,Float64},a)
@test  t == TensorValue{2,2,Float64,3}(11, 12, 21, 22)
t = convert(TensorValue{2,2,Float32},a)
@test  t == TensorValue{2,2,Float32,3}(11, 12, 21, 22)
m = convert(MMatrix{2,2,Float64},t)
@test m == MMatrix{2}(11.0, 12.0, 21.0, 22.0)
u = convert(NTuple{4,Float64},t)
@test u == tuple(11.0, 12.0, 21.0, 22.0)

v = (VectorValue(1.0,2.0),VectorValue(3.0,4.0),VectorValue(5.0,6.0))
t = TensorValues.tensor_from_columns(v...)
@test isa(t,TensorValue{2,3,Float64})
@test t == TensorValue{2,3,Float64}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
t = TensorValues.tensor_from_rows(v...)
@test isa(t,TensorValue{3,2,Float64})
@test t == TensorValue{3,2,Float64}(1.0, 3.0, 5.0, 2.0, 4.0, 6.0)

# Constructors (SymTensorValue)

s = SymTensorValue( (11,21,22) )
@test isa(s,SymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},s) == [11 21;21 22]
@test_throws MethodError convert(SMatrix{1,2,Int},s)

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

a = [11.0 21.0; NaN 22.0]
@test isa(a,AbstractArray{Float64,2})
@test convert(SymTensorValue{2,Float64},a) == SymTensorValue{2,Float64,3}(11.0, 21.0, 22.0)

# Constructors (SymTracelessTensorValue)
q_none = SymTracelessTensorValue{0, Int, 0}()
q = SymTracelessTensorValue()
@test q == q_none
q = SymTracelessTensorValue{0}()
@test q == q_none
q = SymTracelessTensorValue{0}()
@test q == q_none
q = SymTracelessTensorValue(Tuple{}())
@test q == q_none
q = SymTracelessTensorValue{0}(Tuple{}())
@test q == q_none

q_zero = SymTracelessTensorValue{1,Int}(NTuple{0,Int}())
q = SymTracelessTensorValue{1}(Tuple{}())
@test q == q_zero
q = SymTracelessTensorValue{1,Int}(Tuple{}())
@test q == q_zero

q = rand(SymTracelessTensorValue{0,Int})
@test eltype(q) == Int
@test eltype(typeof(q)) == Int

q = SymTracelessTensorValue( (11,21) )
@test isa(q,SymTracelessTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [11 21;21 -11]
@test_throws MethodError convert(SMatrix{1,2,Int},q)

q = SymTracelessTensorValue{2,Int,3}( (11,21) )
@test isa(q,SymTracelessTensorValue{2,Int,3})
@test convert(SMatrix{2,2,Int},q) == [11 21;21 -11]

q = SymTracelessTensorValue{2,Int,3}(11,21)
@test isa(q,SymTracelessTensorValue{2,Int,3})
@test convert(SMatrix{2,2,Int},q) == [11 21;21 -11]

q = SymTracelessTensorValue(11,21)
@test isa(q,SymTracelessTensorValue{2,Int})
@test convert(SMatrix{2,2,Float64},q) == [11.0 21.0;21.0 -11.0]

q = SymTracelessTensorValue{2}( (11,21) )
@test isa(q,SymTracelessTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [11 21;21 -11]

q = SymTracelessTensorValue{2}(11,21)
@test isa(q,SymTracelessTensorValue{2,Int})
@test convert(SMatrix{2,2,Float64},q) == [11.0 21.0;21.0 -11.0]

q = SymTracelessTensorValue{2,Int}( (11,21) )
@test isa(q,SymTracelessTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [11 21;21 -11]

q = SymTracelessTensorValue{2,Float64}(11,21)
@test isa(q,SymTracelessTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},q) == [11.0 21.0;21.0 -11.0]

q = SymTracelessTensorValue{0,Int}( () )
@test isa(q,SymTracelessTensorValue{0,Int})
@test convert(SMatrix{0,0,Int},q) == Array{Any,2}(undef,0,0)

q = SymTracelessTensorValue{0,Int}()
@test isa(q,SymTracelessTensorValue{0,Int})
@test convert(SMatrix{0,0,Int},q) == Array{Any,2}(undef,0,0)

q = SymTracelessTensorValue{1,Int}( () )
@test isa(q,SymTracelessTensorValue{1,Int})
@test convert(SMatrix{1,1,Int},q) == zeros(Int,1,1)

q = SymTracelessTensorValue{1,Int}()
@test isa(q,SymTracelessTensorValue{1,Int})
@test convert(SMatrix{1,1,Int},q) == zeros(Int,1,1)

q = SymTracelessTensorValue(11,21.0)
@test isa(q,SymTracelessTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},q) == [11.0 21.0;21.0 -11.0]

q = SymTracelessTensorValue{2}(11,21.0)
@test isa(q,SymTracelessTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},q) == [11.0 21.0;21.0 -11.0]

q = SymTracelessTensorValue{2,Int}(11,21.0)
@test isa(q,SymTracelessTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [11.0 21.0;21.0 -11.0]

a = [11.0 21.0; NaN NaN]
@test isa(a,AbstractArray{Float64,2})
t = convert(SymTracelessTensorValue{2,Float64},a)
@test  t == SymTracelessTensorValue{2,Float64,3}(11.0, 21.0)
m = convert(MMatrix{2,2,Float64},t)
@test m == MMatrix{2}(11.0, 21.0, 21.0, -11.0)
u = convert(NTuple{3,Float64},t)
@test u == tuple(11.0, 21.0, -11.0)

# Constructors (SkewSymTensorValue)

q_none = SkewSymTensorValue{0, Int, 0}()
q = SkewSymTensorValue()
@test q == q_none
q = SkewSymTensorValue{0}()
@test q == q_none
q = SkewSymTensorValue(Tuple{}())
@test q == q_none
q = SkewSymTensorValue{0}(Tuple{}())
@test q == q_none

q_zero = SkewSymTensorValue{1,Int}(NTuple{0,Int}())
q = SkewSymTensorValue{1}()
@test q == q_zero
q = SkewSymTensorValue{1}(Tuple{}())
@test q == q_zero
q = SkewSymTensorValue{1,Int}(Tuple{}())
@test q == q_zero

q = rand(SkewSymTensorValue{0,Int})
@test eltype(q) == Int
@test eltype(typeof(q)) == Int

q = SkewSymTensorValue( (12,) )
@test isa(q,SkewSymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [0 12;-12 0]
@test_throws MethodError convert(SMatrix{1,2,Int},q)

q = SkewSymTensorValue{2,Int,1}( (12,) )
@test isa(q,SkewSymTensorValue{2,Int,1})
@test convert(SMatrix{2,2,Int},q) == [0 12;-12 0]

q = SkewSymTensorValue{2,Int,1}(12)
@test isa(q,SkewSymTensorValue{2,Int,1})
@test convert(SMatrix{2,2,Int},q) == [0 12;-12 0]

q = SkewSymTensorValue( (12.,13,23) )
@test isa(q,SkewSymTensorValue{3,Float64})
@test convert(SMatrix{3,3,Int},q) == [0 12 13; -12 0 23; -13 -23 0]

q = SkewSymTensorValue(12)
@test isa(q,SkewSymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [0 12;-12 0]

q = SkewSymTensorValue{3}( (12.,13,23) )
@test isa(q,SkewSymTensorValue{3,Float64})
@test convert(SMatrix{3,3,Int},q) == [0 12 13; -12 0 23; -13 -23 0]

q = SkewSymTensorValue{2}( 12 )
@test isa(q,SkewSymTensorValue{2,Int})
@test convert(SMatrix{2,2,Int},q) == [0 12;-12 0]

q = SkewSymTensorValue{3,Int}( (12.,13,23) )
@test isa(q,SkewSymTensorValue{3,Int})
@test convert(SMatrix{3,3,Int},q) == [0 12 13; -12 0 23; -13 -23 0]

q = SkewSymTensorValue{2,Float64}(12)
@test isa(q,SkewSymTensorValue{2,Float64})
@test convert(SMatrix{2,2,Int},q) == [0 12;-12 0]

q = SkewSymTensorValue{3,Int,3}( (12.,13,23) )
@test isa(q,SkewSymTensorValue{3,Int,3})
@test convert(SMatrix{3,3,Int},q) == [0 12 13; -12 0 23; -13 -23 0]

q = SkewSymTensorValue{0,Int}( () )
@test isa(q,SkewSymTensorValue{0,Int})
@test convert(SMatrix{0,0,Int},q) == Array{Any,2}(undef,0,0)

q = SkewSymTensorValue{0,Int}()
@test isa(q,SkewSymTensorValue{0,Int})
@test convert(SMatrix{0,0,Int},q) == Array{Any,2}(undef,0,0)

q = SkewSymTensorValue{1,Int}( () )
@test isa(q,SkewSymTensorValue{1,Int})
@test convert(SMatrix{1,1,Int},q) == zeros(Int,1,1)

q = SkewSymTensorValue{1,Int}()
@test isa(q,SkewSymTensorValue{1,Int})
@test convert(SMatrix{1,1,Int},q) == zeros(Int,1,1)

q = SkewSymTensorValue(12.0)
@test isa(q,SkewSymTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},q) == [0. 12.;-12. 0.]

q = SkewSymTensorValue{2}(12.0)
@test isa(q,SkewSymTensorValue{2,Float64})
@test convert(SMatrix{2,2,Float64},q) == [0. 12.;-12. 0.]

q = SkewSymTensorValue{2,Int}(12.)
@test isa(q,SkewSymTensorValue{2,Int})
@test convert(SMatrix{2,2,Float64},q) == [0 12.;-12. 0]

a = [NaN 12.0; NaN NaN]
@test isa(a,AbstractArray{Float64,2})
t = convert(SkewSymTensorValue{2,Float64},a)
@test  t == SkewSymTensorValue{2,Float64,3}(12.0)
m = convert(MMatrix{2,2,Float64},t)
@test m == MMatrix{2}(0, -12.0, 12.0, 0)
u = convert(NTuple{1,Float64},t)
@test u == tuple(12.0)

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

# Constructors (SymNTensorValue)

# A fully symmetric Nth-order tensor stores the binomial(D+N-1,N) components of
# the non-decreasing indices, in lexicographic order
sn = SymNTensorValue{Tuple{2,2}}((1,2,3))
@test isa(sn,SymNTensorValue{Tuple{2,2},Int,2,3})
sn = SymNTensorValue{Tuple{3,3}}((1,2,3,4,5,6))
@test isa(sn,SymNTensorValue{Tuple{3,3},Int,2,6})
sn = SymNTensorValue{Tuple{2,2,2,2}}((1,2,3,4,5))
@test isa(sn,SymNTensorValue{Tuple{2,2,2,2},Int,4,5})

# Homogeneous NTuple, with and without an explicit component type
@test isa(SymNTensorValue{Tuple{2,2,2}}((1,2,3,4)),
          SymNTensorValue{Tuple{2,2,2},Int,3,4})
@test isa(SymNTensorValue{Tuple{2,2,2},Float64}((1,2,3,4)),
          SymNTensorValue{Tuple{2,2,2},Float64,3,4})
@test isa(SymNTensorValue{Tuple{2,2,2},Float64,3}((1,2,3,4)),
          SymNTensorValue{Tuple{2,2,2},Float64,3,4})
@test isa(SymNTensorValue{Tuple{2,2,2},Float64,3,4}((1.,2.,3.,4.)),
          SymNTensorValue{Tuple{2,2,2},Float64,3,4})

# Heterogeneous Tuple: promoted, or converted to the requested type
@test isa(SymNTensorValue{Tuple{2,2,2}}((1,2.0,3,4)),
          SymNTensorValue{Tuple{2,2,2},Float64,3,4})
@test isa(SymNTensorValue{Tuple{2,2,2},Int}((1,2.0,3,4)),
          SymNTensorValue{Tuple{2,2,2},Int,3,4})

# Vararg
@test SymNTensorValue{Tuple{2,2,2}}(1,2.0,3,4)     == SymNTensorValue{Tuple{2,2,2}}((1.0,2.0,3.0,4.0))
@test SymNTensorValue{Tuple{2,2,2},Int}(1,2.0,3,4) == SymNTensorValue{Tuple{2,2,2}}((1,2,3,4))

# Wrong number of components
@test_throws AssertionError SymNTensorValue{Tuple{3,3}}((1,2,3))
@test_throws AssertionError SymNTensorValue{Tuple{3,3},Int}((1,2,3))
@test_throws AssertionError SymNTensorValue{Tuple{3,3},Int,2}(1,2,3)

# Dimensions that do not all share the same value are not a symmetric shape
@test_throws ErrorException SymNTensorValue{Tuple{2,3}}((1,2,3,4,5,6))
@test_throws ErrorException SymNTensorValue{Tuple{2,2,3}}((1,2,3,4))

# Order 0 holds the single scalar component
sn = SymNTensorValue{Tuple{},Int}((5,))
@test num_indep_components(sn) == 1
@test_throws AssertionError SymNTensorValue{Tuple{},Int}()

# A zero dimension leaves no component at all
sn = SymNTensorValue{Tuple{0,0}}()
@test isa(sn,SymNTensorValue{Tuple{0,0},Int,2,0})
@test num_indep_components(sn) == 0
@test num_components(sn) == 0
@test isa(SymNTensorValue{Tuple{0,0,0},Float64}(), SymNTensorValue{Tuple{0,0,0},Float64,3,0})

# Every fully symmetric tensor is an AbstractSymMultiValue, second-order ones
# included -- symmetry ij<->ji being full symmetry of a pair of indices
@test SymNTensorValue{Tuple{2,2,2}}(1,2,3,4) isa AbstractSymMultiValue
@test SymTensorValue{2}(1,2,3) isa AbstractSymMultiValue
@test SymTensorValue{2}(1,2,3) isa AbstractSymTensorValue
@test AbstractSymTensorValue{2} <: AbstractSymMultiValue
# Symmetry does not fix the stored component count: this one is traceless too
@test SymTracelessTensorValue{2}(1,2) isa AbstractSymMultiValue

# Tensors of a different symmetry class are not
@test !(SkewSymTensorValue{2}(1) isa AbstractSymMultiValue)
# Minor and major symmetry is weaker than full symmetry
@test !(zero(SymFourthOrderTensorValue{2,Int}) isa AbstractSymMultiValue)
@test !(zero(TensorValue{2,2,Int}) isa AbstractSymMultiValue)

# From a full array, which only its non-decreasing indices are read from
for sn in (SymNTensorValue{Tuple{2,2}}(1:3...),
           SymNTensorValue{Tuple{4,4}}(1:10...),
           SymNTensorValue{Tuple{2,2,2,2}}(1:5...),
           SymNTensorValue{Tuple{3,3,3}}(1:10...))
  A = get_array(sn)
  D, N = size(A,1), ndims(A)
  @test SymNTensorValue{NTuple{N,D}}(A) == sn
end

# The full symmetrization of δ⊗δ is fully symmetric, so the round trip through a
# full array reproduces it exactly
kd(i,j) = i == j ? 1 : 0
A_symδδ = [ kd(i,j)*kd(k,l) + kd(i,k)*kd(j,l) + kd(i,l)*kd(j,k) for i in 1:3, j in 1:3, k in 1:3, l in 1:3 ]
symδδ = SymNTensorValue{Tuple{3,3,3,3}}(A_symδδ)
@test get_array(symδδ) == A_symδδ
@test symδδ[1,1,2,2] == symδδ[1,2,1,2] == 1
@test symδδ[1,1,1,1] == 3

# δ⊗δ alone is symmetric in `ij` and in `kl` separately, not under any permutation
A_δδ = [ kd(i,j)*kd(k,l) for i in 1:3, j in 1:3, k in 1:3, l in 1:3 ]
@test A_δδ[1,1,2,2] != A_δδ[1,2,1,2]
@test get_array(SymNTensorValue{Tuple{3,3,3,3}}(A_δδ)) != A_δδ

# Constructors (VectorValue)

a = SVector(1)
g = VectorValue(a)
@test isa(g,VectorValue{1,Int})
@test convert(SVector{1,Int},g) == [1,]
@test_throws MethodError convert(SVector{2,Int},g)
@test convert(SVector{1,Int},g) === SArray(g)
@test convert(MVector{1,Int},g) ==  MArray(g)
@test g === MultiValue(a)

a = SVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test convert(SVector{4,Int},g) == [1,2,3,4]

a = MVector(1,2,3,4)
g = VectorValue(a)
@test isa(g,VectorValue{4,Int})
@test convert(MVector{4,Int},g) == [1,2,3,4]
@test g === MultiValue(a)

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

a = [1.0, 2.0]
@test isa(a,AbstractArray{Float64,1})
t = convert(VectorValue{2,Float64},a)
@test  t == VectorValue{2,Float64}(1.0, 2.0)
m = convert(MVector{2,Float64},t)
@test m == MVector{2}(1.0, 2.0)
u = convert(NTuple{2,Float64},t)
@test u == tuple(1.0, 2.0)

# Initializers

z = zero(TensorValue{3,3,Int,9})
z2= zero(z)
@test z == z2
@test isa(z,TensorValue{3,3,Int,9})
@test convert(SMatrix{3,3,Int},z) == zeros(Int,(3,3))

z = zero(SymTensorValue{3,Int})
z2= zero(z)
@test z == z2
@test isa(z,SymTensorValue{3,Int,6})
@test convert(SMatrix{3,3,Int},z) == zeros(Int,(3,3))

z = zero(SymTracelessTensorValue{3,Int})
z2= zero(z)
@test z == z2
@test isa(z,SymTracelessTensorValue{3,Int,6})
@test convert(SMatrix{3,3,Int},z) == zeros(Int,(3,3))

z = zero(ThirdOrderTensorValue{3,3,3,Int,27})
z2= zero(z)
@test z == z2
@test isa(z,ThirdOrderTensorValue{3,3,3,Int,27})
@test Tuple(z) == Tuple(zeros(Int,(27)))

z = zero(SymFourthOrderTensorValue{2,Int})
z2= zero(z)
@test z == z2
@test isa(z,SymFourthOrderTensorValue{2,Int,9})
@test Tuple(z) == Tuple(zeros(Int,(9)))

z = zero(NTensorValue{Tuple{3,3,3,3},Int,4,81})
z2= zero(z)
@test z == z2
@test isa(z,NTensorValue{Tuple{3,3,3,3},Int,4,81})
@test Tuple(z) == Tuple(zeros(Int,(81)))

z = zero(VectorValue{3,Int})
z2= zero(z)
@test z == z2
@test isa(z,VectorValue{3,Int})
@test convert(SVector{3,Int},z) == zeros(Int,3)

z = one(TensorValue{3,3,Int,9})
@test isa(z,TensorValue{3,3,Int,9})
@test convert(SMatrix{3,3,Int},z) == [1 0 0; 0 1 0; 0 0 1]
s = one(z)
@test convert(SMatrix{3,3,Int},s) == [1 0 0; 0 1 0; 0 0 1]

z = one(TensorValue{2,3,Int,6})
@test isa(z,TensorValue{2,3,Int,6})
@test convert(SMatrix{2,3,Int},z) == [1 0 0; 0 1 0]

z = one(SymFourthOrderTensorValue{2,Int})
@test isa(z,SymFourthOrderTensorValue{2})
@test Tuple(z) == (1.,0.,0., 0.,0.5,0., 0.,0.,1.)

r = rand(VectorValue{3,Int})
@test isa(r,VectorValue{3,Int})
@test r ≠ rand(typeof(r))

r = rand(TensorValue{3,3,Int})
@test isa(r,TensorValue{3,3,Int,9})
@test r ≠ rand(typeof(r))

r = rand(ThirdOrderTensorValue{3,3,3,Int})
@test isa(r,ThirdOrderTensorValue{3,3,3,Int,27})
@test r ≠ rand(typeof(r))

r = rand(SymTensorValue{3,Int})
@test isa(r,SymTensorValue{3,Int,6})
@test r ≠ rand(typeof(r))

r = rand(SymTracelessTensorValue{3,Int})
@test isa(r,SymTracelessTensorValue{3,Int,6})
@test r ≠ rand(typeof(r))

r = rand(SymFourthOrderTensorValue{3,Int})
@test isa(r,SymFourthOrderTensorValue{3,Int,36})
@test r ≠ rand(typeof(r))

r = rand(NTensorValue{Tuple{3,3,3,3},Int})
@test isa(r,NTensorValue{Tuple{3,3,3,3},Int,4,81})
@test r ≠ rand(typeof(r))

sn3 = SymNTensorValue{Tuple{3,3,3}}(1:10...)
z = zero(sn3)
@test z == SymNTensorValue{Tuple{3,3,3}}(ntuple(_->0,10))
@test iszero(z)
@test isa(z,SymNTensorValue{Tuple{3,3,3},Int,3,10})

# Only orders 0 and 2 have a multiplicative neutral element
@test one(SymNTensorValue{Tuple{},Int}) === 1
@test one(SymNTensorValue{Tuple{3,3},Int}) == SymNTensorValue{Tuple{3,3}}(1,0,0,1,0,1)
@test get_array(one(SymNTensorValue{Tuple{3,3},Int})) == [1 0 0; 0 1 0; 0 0 1]
@test_throws ErrorException one(SymNTensorValue{Tuple{3,3,3},Int})
@test_throws ErrorException one(SymNTensorValue{Tuple{3,3,3,3},Int})

# Conversions

a = @SVector ones(Int,3)
b = convert(VectorValue{3,Int},a)
@test isa(b,VectorValue{3,Int})

a = ones(Int,3)
b = convert(VectorValue{3,Int},a)
@test isa(b,VectorValue{3,Int})

ab = Vector(b)
@test ab isa Vector{Int}
@test ab == a
ab = Vector{Float64}(b)
@test ab isa Vector{Float64}
@test ab == a

a = ones(Int,1)
b = convert(VectorValue{1,Int},a)
@test isa(b,VectorValue{1,Int})

a = (1,2,2,1,3,2)
V = TensorValue{3,2,Int,6}
b = convert(V,a)
@test isa(b,V)
b = V[a,a,a,]
@test isa(b,Vector{V})

a = ones(Int,(2,3))
b = convert(TensorValue{2,3,Int},a)
@test isa(b,TensorValue{2,3,Int})

ab = Matrix(b)
@test ab isa Matrix{Int}
@test ab == a
ab = Matrix{Float64}(b)
@test ab isa Matrix{Float64}
@test ab == a

a = (11,21,22)
V = SymTensorValue{2,Int,3}
b = convert(V,a)
@test isa(b,V)
b = V[a,a,a,]
@test isa(b,Vector{V})

a = (11,21)
V = SymTracelessTensorValue{2,Int,3}
b = convert(V,a)
@test isa(b,V)
b = V[a,a,a]
@test isa(b,Vector{V})

a = ones(Int,(2,3,1))
b = convert(ThirdOrderTensorValue{2,3,1,Int},a)
@test isa(b,ThirdOrderTensorValue{2,3,1,Int})

ab = Array(b)
@test ab isa Array{Int}
@test ab == a
ab = Array{Float64}(b)
@test ab isa Array{Float64}
@test ab == a

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

a = ones(Int,(2,3,1,5))
b = convert(NTensorValue{Tuple{2,3,1,5},Int},a)
@test isa(b,NTensorValue{Tuple{2,3,1,5},Int})

# Misc operations on the type itself

V = VectorValue{3,Int}
@test length(V) == 3
@test size(V) == (3,)

V = VectorValue{3}
@test length(V) == 3
@test size(V) == (3,)

# Too abstract types error instead of silently giving a wrong answer

@test_throws "has no element type" eltype(VectorValue{3})
@test_throws "the size information is needed" length(VectorValue)
@test_throws "has no size information" size(VectorValue)
@test_throws "the size and eltype parameters" zero(VectorValue{3})
@test_throws "parametric argument) is needed" num_components(VectorValue)

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
@test convert(SArray{Tuple{2,2,2},Int},t) === SArray(t)
@test convert(MArray{Tuple{2,2,2},Int},t) == MArray(t)
@test t === MultiValue(a)

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

# Fourth order tensor

a = SArray{Tuple{2,2,2,2}}(1:16...)
t = NTensorValue(a)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Int})
@test convert(SArray{Tuple{2,2,2,2},Int},t) === SArray(t)
@test convert(MArray{Tuple{2,2,2,2},Int},t) == MArray(t)
@test t === MultiValue(a)

t = NTensorValue(Val(4),1:16...)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Int})

t = NTensorValue{Tuple{2,2,2,2}}(1:16...)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Int})

t = NTensorValue{Tuple{2,2,2,2},Float64}(1:16...)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Float64})

t = NTensorValue(Val(4),1,2.0,3:16...)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Float64})

t = NTensorValue{Tuple{2,2,2,2}}(1,2.0,3:16...)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Float64})

t = NTensorValue{Tuple{2,2,2,2},Int}(1.,2.0,3:16...)
@test isa(t,NTensorValue{Tuple{2,2,2,2},Int})


# mutable, Mutable, implicit convert

v = VectorValue(3,2,1)
m = mutable(v)
@test m == get_array(v) == SArray(v)
@test isa(m,MVector)

v = TensorValue{2,3}(1,2,3,4,5,6)
m = mutable(v)
@test m == get_array(v) == SArray(v)
@test isa(m,MMatrix)

v = SymTensorValue{2}(1,2,3)
m = mutable(v)
@test m == get_array(v) == SArray(v)
@test isa(m,MMatrix)

v = SymTracelessTensorValue{2}(1,2)
m = mutable(v)
@test m == get_array(v) == SArray(v)
@test isa(m,MMatrix)

v = ThirdOrderTensorValue{2,1,3}(1:6...)
m = mutable(v)
@test m == get_array(v) == SArray(v)
@test isa(m,MArray)

v = SymFourthOrderTensorValue{2}(1:9...)
m = mutable(v)
@test m == get_array(v)
@test isa(m,MArray)

v = NTensorValue{Tuple{2,1,3,1}}(1:6...)
m = mutable(v)
@test m == get_array(v) == SArray(v)
@test isa(m,MArray)

V = VectorValue{3,Int}
M = Mutable(V)
@test M == MVector{3,Int}
m = rand(M)
v = VectorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ] # convert test
@test M[ v ] == M[ m ] # convert test

V = TensorValue{3,3,Int}
M = Mutable(V)
@test M == MMatrix{3,3,Int,9}
m = rand(M)
v = TensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test M[ v ] == M[ m ]
@test V[ v ] == V[ m ]

V = SymTensorValue{3,Int}
M = Mutable(V)
@test M == MMatrix{3,3,Int,9}
m = rand(M)
v = SymTensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ]

V = SymTracelessTensorValue{3,Int}
M = Mutable(V)
@test M == MMatrix{3,3,Int,9}
m = rand(M)
v = SymTracelessTensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ]

V = SkewSymTensorValue{3,Int}
M = Mutable(V)
@test M == MMatrix{3,3,Int,9}
m = rand(M)
v = SkewSymTensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ]

V = ThirdOrderTensorValue{3,1,2,Int}
M = Mutable(V)
@test M == MArray{Tuple{3,1,2},Int,3,6}
m = rand(M)
v = ThirdOrderTensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ]
@test M[ v ] == M[ m ]

V = SymFourthOrderTensorValue{2,Int}
M = Mutable(V)
@test M == MArray{Tuple{2,2,2,2},Int,4,16}
m = zero(M)
v = SymFourthOrderTensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ]

V = NTensorValue{Tuple{3,1,2,1},Int}
M = Mutable(V)
@test M == MArray{Tuple{3,1,2,1},Int,4,6}
m = rand(M)
v = NTensorValue(m)
@test isa(v,V)
M2 = Mutable(v)
@test M == M2
@test V[ v ] == V[ m ]
@test M[ v ] == M[ m ]

@test_throws ErrorException Mutable(MultiValue) # @abstractmethod

sn3 = SymNTensorValue{Tuple{3,3,3}}(1:10...)
@test isa(get_array(sn3), SArray{Tuple{3,3,3},Int,3,27})
@test isa(SArray(sn3), SArray{Tuple{3,3,3},Int,3,27})
@test isa(Array(sn3), Array{Int,3})
@test Array(sn3) == get_array(sn3)
@test isa(mutable(sn3), MArray{Tuple{3,3,3},Int,3,27})
@test Mutable(sn3) == MArray{Tuple{3,3,3},Int,3,27}
@test convert(SymNTensorValue{Tuple{3,3,3},Float64}, sn3) ==
      SymNTensorValue{Tuple{3,3,3},Float64}(1:10...)
@test convert(SymNTensorValue{Tuple{3,3,3},Int}, sn3) === sn3

# conversion with scalar for Tensors with exactly one physical component
v  = VectorValue(1)
t2 = TensorValue(1)
s  = SymTensorValue(1)
t3 = ThirdOrderTensorValue(1)
f  = SymFourthOrderTensorValue(1)
h5 = NTensorValue(Val(5),1)
@test Float64[ v ]  == [ 1. ]
@test Float64[ t2 ] == [ 1. ]
@test Float64[ s ]  == [ 1. ]
@test Float64[ t3 ] == [ 1. ]
@test Float64[ f ]  == [ 1. ]
@test Float64[ h5 ]  == [ 1. ]
@test VectorValue{1}[ 1 ] == [ v ]
@test TensorValue{1,1,Int,1}[ 1 ] == [ t2 ]
@test SymTensorValue{1,Int}[ 1 ] == [ s ]
@test ThirdOrderTensorValue{1,1,1,Int,1}[ 1 ] == [ t3 ]
@test SymFourthOrderTensorValue{1,Int}[ 1 ] == [ f ]
@test NTensorValue{Tuple{1,1,1,1,1},Int,5,1}[ 1 ] == [ h5 ]

v  = VectorValue(1,2)
st = SymTracelessTensorValue{1}()
sk = SkewSymTensorValue(1)
@test Float64[ st ] == [ 0. ]
@test_throws ErrorException Float64[ v ]
@test_throws ErrorException Float64[ sk ]
@test_throws ErrorException VectorValue{2,Int,2}[ 1 ]
@test_throws ErrorException SymTracelessTensorValue{1,Int,0}[ 1 ]
@test_throws ErrorException SkewSymTensorValue{1,Int,1}[ 1 ]
@test_throws ErrorException NTensorValue{Tuple{1,2,1,1},Int,2}[ 1 ]

# component number APIs
@test num_components(Int) == 1
@test num_components(Float64) == 1
@test num_components(1.0) == 1
@test num_components(1) == 1
@test num_components(VectorValue{3,Float64}) == 3
@test num_components(AbstractSymTensorValue{2}) == 4
@test num_components(VectorValue(1,2,3)) == 3
@test num_components(TensorValue(1,2,3,4)) == 4
@test num_components(SymTensorValue(1,2,3)) == 4
@test num_components(SymTracelessTensorValue(1,2)) == 4
@test num_components(ThirdOrderTensorValue(111,112,121,122,211,212,221,222)) == 8
@test num_components(SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)) == 16
@test num_components(NTensorValue{Tuple{2,3,1,5}}) == 2*3*1*5
@test num_components(SymNTensorValue{Tuple{3,3,3},Int}) == 27

@test num_indep_components(Int) == 1
@test num_indep_components(Float64) == 1
@test num_indep_components(1.0) == 1
@test num_indep_components(1) == 1
@test num_indep_components(VectorValue{0}) == 0
@test num_indep_components(VectorValue{3}) == 3
@test num_indep_components(VectorValue(1,2,3)) == 3
@test num_indep_components(TensorValue{2,2}) == 4
@test num_indep_components(TensorValue(1,2,3,4)) == 4
@test num_indep_components(SymTensorValue{0}) == 0
@test num_indep_components(SymTensorValue{2}) == 3
@test num_indep_components(SymTensorValue(1,2,3)) == 3
@test num_indep_components(SymTracelessTensorValue{0}) == 0
@test num_indep_components(SymTracelessTensorValue{1}) == 0
@test num_indep_components(SymTracelessTensorValue{2}) == 2
@test num_indep_components(SymTracelessTensorValue(1,2)) == 2
@test num_indep_components(ThirdOrderTensorValue{2,2,2}) == 8
@test num_indep_components(ThirdOrderTensorValue(111,112,121,122,211,212,221,222)) == 8
@test num_indep_components(SymFourthOrderTensorValue{2}) == 9
@test num_indep_components(SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)) == 9
@test num_indep_components(NTensorValue{Tuple{2,3,1,5}}) == 2*3*1*5
@test num_indep_components(SymNTensorValue{Tuple{3,3,3}}) == 10
@test num_indep_components(SymNTensorValue{Tuple{2,2,2,2}}) == 5

@test_throws ErrorException num_components(VectorValue)
@test_throws ErrorException num_components(TensorValue)
@test_throws ErrorException num_components(TensorValue{2})
@test_throws ErrorException num_components(AbstractSymTensorValue)
@test_throws ErrorException num_components(SymTensorValue)
@test_throws ErrorException num_components(SymTracelessTensorValue)
@test_throws ErrorException num_components(ThirdOrderTensorValue{2,2})
@test_throws ErrorException num_components(SymFourthOrderTensorValue)
@test_throws ErrorException num_indep_components(VectorValue)
@test_throws ErrorException num_indep_components(TensorValue)
@test_throws ErrorException num_indep_components(TensorValue{2})
@test_throws ErrorException num_indep_components(AbstractSymTensorValue{2})
@test_throws ErrorException num_indep_components(SymTensorValue)
@test_throws ErrorException num_indep_components(SymTracelessTensorValue)
@test_throws ErrorException num_indep_components(ThirdOrderTensorValue{2,2})
@test_throws ErrorException num_indep_components(SymFourthOrderTensorValue)
@test_throws ErrorException num_indep_components(NTensorValue)

# To vizualise component names:
# tt = TensorValue(1:4...)
# for ci in CartesianIndices(tt)
#   di = tt[ci] # di = data_index(tt,Tuple(ci)...)
#   println(di, " ", indep_components_names(tt)[di]," ", ci)
# end
@test indep_components_names(VectorValue{3}) == ["X","Y","Z"]
@test indep_components_names(VectorValue{4}) == ["1","2","3","4"]
@test indep_components_names(TensorValue{2,2}) == ["XX","YX","XY","YY"]
@test indep_components_names(TensorValue{2,4}) == ["11","21","12","22","13","23","14","24"]
@test indep_components_names(SymTensorValue{2}) == ["XX","XY","YY"]
@test indep_components_names(SymTensorValue{4}) == ["11","12","13","14","22","23","24","33","34","44"]
@test indep_components_names(SymTracelessTensorValue{2}) == ["XX","XY","YY"]
@test indep_components_names(SymTracelessTensorValue{4}) == ["11","12","13","14","22","23","24","33","34","44"]
@test indep_components_names(SkewSymTensorValue{2}) == ["XY"]
@test indep_components_names(SkewSymTensorValue{3}) == ["XY", "XZ", "YZ"]
@test indep_components_names(ThirdOrderTensorValue{2,2,1}) == ["XXX","YXX","XYX","YYX"]
@test indep_components_names(ThirdOrderTensorValue{1,4,1}) == ["111","121","131","141"]
@test indep_components_names(SymFourthOrderTensorValue{2}) == [
 "XXXX", "XXXY", "XXYY", "XYXX", "XYXY", "XYYY", "YYXX", "YYXY", "YYYY"
]
@test indep_components_names(NTensorValue{Tuple{1,2,2,1}}) == ["XXXX","XYXX","XXYX","XYYX"]
@test indep_components_names(NTensorValue{Tuple{1,1,4,1}}) == ["1111","1121","1131","1141"]

@test indep_components_names(SymNTensorValue{Tuple{2,2,2},Int}) == ["XXX","XXY","XYY","YYY"]
@test indep_components_names(SymNTensorValue{Tuple{3,3},Int}) == ["XX","XY","XZ","YY","YZ","ZZ"]
@test indep_components_names(SymNTensorValue{Tuple{4,4},Int}) ==
      indep_components_names(SymTensorValue{4,Int})
# On the fully parameterized type too, which is what the VTK export passes
@test indep_components_names(SymNTensorValue{Tuple{3,3,3},Int,3,10}) ==
      [ join(("X","Y","Z")[i] for i in I) for I in sorted_multisets(3,3) ]

V = MultiValue{Tuple{3,4},Int,2,5}
warn_log = "indep_components_names has not been implemented"
@test (@test_warn  warn_log indep_components_names(V)) == ["1","2","3","4","5"]

a = VectorValue(1,2,3,4)
@test change_eltype(a,Float64) == VectorValue{4,Float64}

a = TensorValue(1,2,3,4)
@test change_eltype(a,Float64) == TensorValue{2,2,Float64,4}
@test change_eltype(TensorValue{2,2,Float64},Int) == TensorValue{2,2,Int}

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
@test change_eltype(SymTensorValue{2,Float64},Int) == SymTensorValue{2,Int}
@test isa(Tuple(a),Tuple)
@test Tuple(a) == a.data
b = Matrix{Int}(undef,2,2)
b[1,1] = a[1,1]
b[1,2] = a[1,2]
b[2,1] = a[2,1]
b[2,2] = a[2,2]
a = SymTensorValue(11,21,22)
bt = SymTensorValue{2,Int}(b)
@test all(bt .== a)

a = SymTracelessTensorValue(11,21)
@test change_eltype(a,Float64) == SymTracelessTensorValue{2,Float64,3}
@test change_eltype(SymTracelessTensorValue{2,Float64},Int) == SymTracelessTensorValue{2,Int}
@test isa(Tuple(a),Tuple)
@test Tuple(a) == a.data
b = Matrix{Int}(undef,2,2)
b[1,1] = a[1,1]
b[1,2] = a[1,2]
b[2,1] = a[2,1]
b[2,2] = a[2,2]
bt = SymTracelessTensorValue{2,Int}(b)
@test all(bt .== a)

a = SymFourthOrderTensorValue(1111,1121,1122, 2111,2121,2122, 2211,2221,2222)
@test change_eltype(a,Float64) == SymFourthOrderTensorValue{2,Float64,9}
@test change_eltype(SymFourthOrderTensorValue{2,Float64},Int) == SymFourthOrderTensorValue{2,Int}
@test isa(Tuple(a),Tuple)
@test Tuple(a) == a.data

sn3 = SymNTensorValue{Tuple{3,3,3}}(1:10...)
@test change_eltype(sn3,Float64) == SymNTensorValue{Tuple{3,3,3},Float64,3,10}
@test change_eltype(SymNTensorValue{Tuple{3,3,3},Int,3,10},Float64) ==
      SymNTensorValue{Tuple{3,3,3},Float64,3,10}
@test isa(Tuple(sn3),Tuple)
@test Tuple(sn3) == sn3.data
@test promote_type(SymNTensorValue{Tuple{3,3,3},Int,3,10},
                   SymNTensorValue{Tuple{3,3,3},Float64,3,10}) ==
      SymNTensorValue{Tuple{3,3,3},Float64,3}


for V in (VectorValue, SymTensorValue, SkewSymTensorValue,
            SymTracelessTensorValue, SymFourthOrderTensorValue)

  VD = V{3,Float32}
  u = rand(VD)
  comp_basis = component_basis(VD)
  dual_basis_rpzs = representatives_of_componentbasis_dual(VD)
  @test comp_basis == component_basis(u)
  @test dual_basis_rpzs == representatives_of_componentbasis_dual(u)

  @test length(comp_basis) == length(dual_basis_rpzs) == num_indep_components(VD)
  @test eltype(comp_basis) <: VD
  @test eltype(dual_basis_rpzs) <: VD
  @test u ≈ sum( indep_comp_getindex(u,i) * Vi for (i,Vi) in enumerate(comp_basis) )
  @test u ≈ V( (u ⊙ Vi for Vi in dual_basis_rpzs)... )
end

end # module TypesTests
