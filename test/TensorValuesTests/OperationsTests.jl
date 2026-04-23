module OperationsTests

using Test: Error
using Test
using Gridap.TensorValues
using Gridap.Arrays
using LinearAlgebra
using Gridap.TensorValues: _eltype
using StaticArrays

# Comparison
a0 = 0
a1 = 1
b = VectorValue(1,3,3)

@test (a0 < b) == true
@test (a0 <= b) == true
@test (a0 == b) == false
@test (a1 < b) == false
@test (a1 <= b) == true


a = VectorValue(1,2,3)
b = VectorValue(1,3,3)

@test (a < b) == true
@test (a <= b) == true
@test (a == b) == false
@test (a >= b) == false
@test (a > b) == false

@test (a < a) == false
@test (a <= a) == true
@test (a == a) == true
@test (a >= a) == true
@test (a > a) == false

@test VectorValue(1,2,3) == VectorValue(1.0,2.0,3.0)
@test VectorValue(1,2,3) == VectorValue(1+0im, 2+0im, 3+0im)
@test VectorValue(1,2,3) ≠ VectorValue(1,2)
@test VectorValue(1,2,3) ≠ SymTensorValue(1,2,3)
@test iszero(VectorValue(1,2,3) - VectorValue(1.0,2.0,3.0))
@test iszero(zero(VectorValue(1,2,3)))
@test isapprox(VectorValue(1,2,3), VectorValue(1.0,2.0,3.0))
@test isapprox(VectorValue(1,2,3), VectorValue(1.0,2.0,3.0), atol=eps(1.0))
@test isapprox(VectorValue(1,2,3), VectorValue(1.0,2.0,3.0), rtol=eps(1.0))
@test !isapprox(VectorValue(1,2,3), VectorValue(1.0,2.0,4.0))
@test_throws DimensionMismatch isapprox(VectorValue(1,2,3), VectorValue(1.0,2.0,3.0,4.0))

a = VectorValue(1,2,3)
b = VectorValue(2,1,6)

@test a==a
@test a ≈ a
@test a!=b
@test [a,a] == [a,a]
@test [a,a] ≈ [a,a]

c = TensorValue(1:9...)
@test !(a==c)
@test_throws DimensionMismatch !(a≈c)
@test !([a,a]==[c,c])
@test !([a,a,a]≈[c,c])
@test_throws DimensionMismatch !([a,a]≈[c,c])
@test_throws ErrorException (a < c)
@test_throws ErrorException (a <= c)

# Addition / subtraction

c = +a
r = VectorValue(1,2,3)
@test c == r

c = -a
r = VectorValue(-1,-2,-3)
@test c == r

c = a + b
r = VectorValue(3,3,9)
@test c == r

c = a - b
r = VectorValue(-1,1,-3)
@test c == r

a = VectorValue{0,Float64}()
b = VectorValue{0,Int}()
@test a == b
@test a ≈ b

c = +a
r = VectorValue{0,Float64}()
@test c == r

c = -a
r = VectorValue{0,Float64}()
@test c == r

c = a + b
r = VectorValue{0,Float64}()
@test c == r

c = a - b
r = VectorValue{0,Float64}()
@test c == r

c = a*2
r = VectorValue{0,Float64}()
@test c == r

c = 2*a
r = VectorValue{0,Float64}()
@test c == r

a = TensorValue(1,2,3,4)
b = TensorValue(5,6,7,8)

c = +a
r = a
@test c==r

c = -a
r = TensorValue(-1,-2,-3,-4)
@test c==r

c = a - b
r = TensorValue(-4, -4, -4, -4)
@test c==r

b  = TensorValue(5.,6,7,8) # promote_rule test
c = a - b
r = TensorValue(-4.,-4,-4,-4)
@test c===r

a = SymTensorValue(1,2,3)
b = SymTensorValue(5,6,7)

c = -a
r = SymTensorValue(-1,-2,-3)
@test c==r

c = a + b
r = SymTensorValue(6,8,10)
@test c==r

c = b - a
r = SymTensorValue(4,4,4)
@test c==r

b = SymTensorValue(5.,6,7)
c = a + b
r = SymTensorValue(6.,8,10)
@test c===r

a = TensorValue(1,2,3,4)
b = SymTensorValue(5,6,7)

c = a + b
d = b + a
r = TensorValue(6,8,9,11)
@test c==r
@test d==r

c = a - b
r = TensorValue(-4,-4,-3,-3)
@test c==r

c = b - a
r = TensorValue(4,4,3,3)
@test c==r

a = SymTracelessTensorValue(1,2)
b = SymTracelessTensorValue(5,6)

c = -a
r = SymTracelessTensorValue(-1,-2)
@test c==r

c = a + b
r = SymTracelessTensorValue(6,8)
@test c==r

c = a - b
r = SymTracelessTensorValue(-4,-4)
@test c==r

b = SymTracelessTensorValue(5.,6)
c = a + b
r = SymTracelessTensorValue(6.,8)
@test c===r

a = SymTensorValue(1,2,3)
b = SymTracelessTensorValue(5,6)

c = a + b
d = b + a
r = SymTensorValue(6,8,-2)
@test c==r
@test d==r

c = a - b
r = SymTensorValue(-4,-4,8)
@test c==r

c = b - a
r = SymTensorValue(4,4,-8)
@test c==r

a = SkewSymTensorValue(1,2,3)
b = SkewSymTensorValue(4,5,6)

c = -a
r = SkewSymTensorValue(-1,-2,-3)
@test c==r

c = a + b
r = SkewSymTensorValue(5,7,9)
@test c==r

c = a - b
r = SkewSymTensorValue(-3,-3,-3)
@test c==r

b = SkewSymTensorValue(4.,5,6)
c = a + b
r = SkewSymTensorValue(5.,7,9)
@test c===r

a = SkewSymTensorValue(1,2,3)
b = SymTracelessTensorValue(1:5...)

c = a + b
d = b + a
r = TensorValue(1, 1, 1, 3, 4, 2, 5, 8, -5)
@test c==r
@test d==r

c = a - b
r = TensorValue(-1,-3,-5,-1,-4,-8,-1,-2,5)
@test c==r

c = b - a
r = TensorValue(1,3,5,1,4,8,1,2,-5)
@test c==r

a = SkewSymTensorValue(1,2,3)
b = SymTensorValue(1:6...)

c = a + b
d = b + a
r = TensorValue(1, 1, 1, 3, 4, 2, 5, 8, 6)
@test c==r
@test d==r

c = a - b
r = TensorValue(-1,-3,-5,-1,-4,-8,-1,-2,-6)
@test c==r

c = b - a
r = TensorValue(1,3,5,1,4,8,1,2,6)
@test c==r


a = SymTensorValue(5,6,7)
b = TensorValue(1,2,3,4)

c = a + b
d = b + a
r = TensorValue(6,8,9,11)
@test c==r
@test d==r

c = a - b
r = TensorValue(4,4,3,3)
@test c==r

c = b - a
r = TensorValue(-4,-4,-3,-3)
@test c==r

a = SymTracelessTensorValue(5,6)
b = TensorValue(1,2,3,4)

c = a + b
d = b + a
r = TensorValue(6,8,9,-1)
@test c==r
@test d==r

c = a - b
r = TensorValue(4,4,3,-9)
@test c==r

c = b - a
r = TensorValue(-4,-4,-3,9)
@test c==r

a = SkewSymTensorValue(1,2,3)
b = TensorValue(1:9...)

c = a + b
d = b + a
r = TensorValue(1, 1, 1, 5, 5, 3, 9, 11, 9)
@test c==r
@test d==r

c = a - b
r = TensorValue(-1,-3,-5,-3,-5,-9,-5,-5,-9)
@test c==r

c = b - a
r = TensorValue(1,3,5,3,5,9,5,5,9)
@test c==r

a = ThirdOrderTensorValue(1:8...)
b = ThirdOrderTensorValue(9:16...)

c = -a
r=ThirdOrderTensorValue(-1,-2,-3,-4,-5,-6,-7,-8)
@test c==r

c = a + b
r=ThirdOrderTensorValue(10,12,14,16,18,20,22,24)
@test c==r

c = a - b
r=ThirdOrderTensorValue(-8,-8,-8,-8,-8,-8,-8,-8)
@test c==r

b = ThirdOrderTensorValue(9.,10:16...)
c = a + b
r = ThirdOrderTensorValue(10.,12,14,16,18,20,22,24)
@test c===r


v = VectorValue(1,2)
t = TensorValue(1,2,3,4)
s = SymTensorValue(1,2,3)
q = SymTracelessTensorValue(1,2)
k = SkewSymTensorValue(1)
r = ThirdOrderTensorValue(1:8...)
f = SymFourthOrderTensorValue(1:9...)

@test_throws ErrorException v+t
@test_throws ErrorException r+v
@test_throws ErrorException r-f
@test_throws ErrorException f-v
@test_throws ErrorException v-s
@test_throws ErrorException q+v
@test_throws ErrorException v+q
@test_throws ErrorException v-q
@test_throws ErrorException q-v
@test_throws ErrorException v-k
@test_throws ErrorException k-v

@test_throws ErrorException q+0
@test_throws ErrorException 0+q
@test_throws ErrorException 0-q
@test_throws ErrorException q-0

# Multiplication / division

@test_throws ErrorException v*t
@test_throws ErrorException r*v
@test_throws ErrorException r/f
@test_throws ErrorException f/v
@test_throws ErrorException v/s
@test_throws ErrorException q*v
@test_throws ErrorException v*q
@test_throws ErrorException q*s
@test_throws ErrorException s*q
@test_throws ErrorException s*s
@test_throws ErrorException q*q
@test_throws ErrorException v/q
@test_throws ErrorException q/v
@test_throws ErrorException v*k
@test_throws ErrorException k*v
@test_throws ErrorException v/k
@test_throws ErrorException k/v

# Matrix Division

a = VectorValue(1,2,3)

t = one(TensorValue{3,3,Int})
c = t\a
@test c == a

st = one(SymTensorValue{3,Int})
c = st\a
@test c == a

@test_throws ErrorException one(SkewSymTensorValue{3,Int})
@test_throws ErrorException one(MultiValue)

# Operations by a scalar

t = TensorValue(1,2,3,4,5,6,7,8,9)
st = SymTensorValue(1,2,3,5,6,9)
qt = SymTracelessTensorValue(1,2,3,5,6)
sk = SkewSymTensorValue(1,2,3)
s4ot = one(SymFourthOrderTensorValue{2,Int})
a = VectorValue(1,2,3)

c = 2 * a
@test isa(c,VectorValue{3,Int})
r = VectorValue(2,4,6)
@test c == r

c = a * 2
@test isa(c,VectorValue{3,Int})
r = VectorValue(2,4,6)
@test c == r

c = 2 + a
@test isa(c,VectorValue{3,Int})
r = VectorValue(3,4,5)
@test c == r

c = a + 2
@test isa(c,VectorValue{3,Int})
@test c == r

c = a / 2
@test isa(c,VectorValue{3,Float64})
r = VectorValue(1/2,1.0,3/2)
@test c == r

c = 2 * t
@test isa(c,TensorValue{3})
r = TensorValue(2, 4, 6, 8, 10, 12, 14, 16, 18)
@test c == r

c = t * 2
@test isa(c,TensorValue{3})
r = TensorValue(2, 4, 6, 8, 10, 12, 14, 16, 18)
@test c == r

c = t / 2
@test isa(c,TensorValue{3})
r = TensorValue(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)
@test c == r

c = t + 2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(3, 4, 5, 6, 7, 8, 9, 10, 11)
@test c == r


c = 2 * st
@test isa(c,SymTensorValue{3})
r = SymTensorValue(2,4,6,10,12,18)
@test c == r

c = st * 2
@test isa(c,SymTensorValue{3})
r = SymTensorValue(2,4,6,10,12,18)
@test c == r

c = st / 2
@test isa(c,SymTensorValue{3})
r = SymTensorValue(.5,1,1.5,2.5,3,4.5)
@test c == r

c = st + 2
@test isa(c,SymTensorValue{3})
r = SymTensorValue(3,4,5,7,8,11)
@test c == r


c = 2 * qt
@test isa(c,SymTracelessTensorValue{3})
r = SymTracelessTensorValue(2,4,6,10,12)
@test c == r

c = qt * 2
@test isa(c,SymTracelessTensorValue{3})
r = SymTracelessTensorValue(2,4,6,10,12)
@test c == r

c = qt / 2
@test isa(c,SymTracelessTensorValue{3})
r = SymTracelessTensorValue(.5,1,1.5,2.5,3)
@test c == r

c = 2 * sk
@test isa(c,SkewSymTensorValue{3})
r = SkewSymTensorValue(2,4,6)
@test c == r

c = sk * 2
@test isa(c,SkewSymTensorValue{3})
r = SkewSymTensorValue(2,4,6)
@test c == r

c = sk / 2
@test isa(c,SkewSymTensorValue{3})
r = SkewSymTensorValue(.5,1,1.5)
@test c == r

c = 2 * s4ot
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(2,0,0, 0,1,0, 0,0,2)
@test c == r

c = s4ot * 2
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(2,0,0, 0,1,0, 0,0,2)
@test c == r

c = s4ot / 2
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(.5,0,0, 0,.25,0, 0,0,.5)
@test c == r

c = s4ot + 0
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(1,0,0, 0,.5,0, 0,0,1)
@test c == r

# Dot product (simple contraction)

a = VectorValue(1,2,3)
b = VectorValue(2,1,6)

t = TensorValue(1,2,3,4,5,6,7,8,9)
s = TensorValue(9,8,3,4,5,6,7,2,1)
st  = SymTensorValue(1,2,3,5,6,9)
st2 = SymTensorValue(9,6,5,3,2,1)
qt  = SymTracelessTensorValue(1,2,3,5,6)
qt2 = SymTracelessTensorValue(9,6,5,3,2)
sk  = SkewSymTensorValue(1,2,3)
sk2 = SkewSymTensorValue(3,2,1)

c = a ⋅ b
@test isa(c,Int)
@test c == 2+2+18

c = t ⋅ a
@test isa(c,VectorValue{3,Int})
r = VectorValue(30,36,42)
@test c == r

c = st ⋅ a
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,42)
@test c == r

c = qt ⋅ a
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,-3)
@test c == r

c = sk ⋅ a
@test isa(c,VectorValue{3,Int})
r = VectorValue(8,8,-8)
@test c == r

c = s ⋅ t
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(38,24,18,98,69,48,158,114,78)
@test c == r

c = st ⋅ st2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 78, 108, 18, 39, 54, 12, 26, 36)
@test c == r

c = qt ⋅ qt2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 78, 33, 18, 39, 24, -27, -52, 99)
@test c == r

c = sk ⋅ sk2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(-7, -6, 9, -2, -6, -6, 1, -2, -7)
@test c == r

c = st ⋅ qt2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 78, 108, 18, 39, 54, -27, -52, -81)
@test c == r

c = qt2 ⋅ st
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 18, -27, 78, 39, -52, 108, 54, -81)
@test c == r

c = st ⋅ sk2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(-12, -27, -36, 0, 0, 0, 4, 9, 12)
@test c == r

c = sk2 ⋅ st
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(12, 0, -4, 27, 0, -9, 36, 0, -12)
@test c == r

c = a ⋅ st
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,42)
@test c == r

c = a ⋅ qt
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,-3)
@test c == r

c = a ⋅ sk
@test isa(c,VectorValue{3,Int})
r = VectorValue(-8,-8,8)
@test c == r

a1 = VectorValue(1,0)
b1 = VectorValue(1,2)

t1 = ThirdOrderTensorValue{2,2,1}(1,2,3,4)
t2 = TensorValue(1,0,0,1)
t3 = TensorValue(1,2,0,0)
t4 = ThirdOrderTensorValue{2,2,2}(1,2,3,4,5,6,7,8)
t5 = TensorValue{2,1}(1,2)

c = a1 ⋅ t1
@test isa(c,TensorValue{2,1,Int})
r = TensorValue{2,1}(1,3)
@test c == r

c = b1 ⋅ t1
@test isa(c,TensorValue{2,1,Int})
r = TensorValue{2,1}(5,11)
@test c == r

c = t2 ⋅ t1
@test isa(c,ThirdOrderTensorValue{2,2,1,Int,4})
r = ThirdOrderTensorValue{2,2,1}(1,2,3,4)
@test c == r

c = t3 ⋅ t1
@test isa(c,ThirdOrderTensorValue{2,2,1,Int,4})
r = ThirdOrderTensorValue{2,2,1}(1,2,3,6)
@test c == r

c = t4 ⋅ t3
@test isa(c,ThirdOrderTensorValue{2,2,2,Int,8})
r = ThirdOrderTensorValue{2,2,2}(11,14,17,20,0,0,0,0)
@test c == r

c = t4 ⋅ t5
@test isa(c,ThirdOrderTensorValue{2,2,1,Int,4})
r = ThirdOrderTensorValue{2,2,1}(11,14,17,20)
@test c == r

x = VectorValue{0,Float64}()
G = TensorValue{0,2,Float64,0}()
@test x⋅G == VectorValue(0,0)

# Inner product (full contraction)

c = 2 ⊙ 3
@test c == 6

c = a ⊙ b
@test isa(c,Int)
@test c == 2+2+18

c = inner(t,s)
@test isa(c,Int)
@test c == 185

c = inner(st,st2)
@test c == st ⊙ st2
@test isa(c,Int)
@test c == inner(TensorValue(get_array(st)),TensorValue(get_array(st2)))

c = inner(qt,qt2)
@test c == qt ⊙ qt2
@test isa(c,Int)
@test c == inner(TensorValue(get_array(qt)),TensorValue(get_array(qt2)))

t1 = TensorValue{2,3}(1:6...)
t2 = TensorValue{2,3}(10:15...)
@test inner(t1,t1) == 91
@test double_contraction(t1,t2) == inner(t1,t2)

c = inner(t4,t4)
@test c == t4 ⊙ t4
@test isa(c,Int)
@test c == inner(t4[:],t4[:])

@test_throws ErrorException inner(a,t)
@test_throws ErrorException inner(s,a)
@test_throws ErrorException inner(a,q)
@test_throws ErrorException inner(a,t4)

# Reductions

a = VectorValue(1,2,3)

c = sum(a)
@test isa(c,Int)
@test c == 1+2+3

c = maximum(a)
@test isa(c,Int)
@test c == 3

c = minimum(a)
@test isa(c,Int)
@test c == 1

# Outer product (aka dyadic product)

a = VectorValue(1,2,3)
e = VectorValue(2,5)

c = outer(2,3)
c = 2 ⊗ 3
@test c == 6

r = VectorValue(2,4,6)
c = outer(2,a)
c = 2 ⊗ a
@test isa(c,VectorValue{3,Int})
@test c == r

c = outer(a,2)
@test isa(c,VectorValue{3,Int})
@test c == r

c = outer(a,e)
c = a ⊗ e
@test isa(c,TensorValue{3,2,Int})
r = TensorValue{3,2,Int}(2,4,6,5,10,15)
@test c == r

e = VectorValue(10,20,30)
k = TensorValue{4,1}(1,2,3,4)
c = outer(e,k)
@test c == ThirdOrderTensorValue{3,4,1}(10, 20, 30, 20, 40, 60, 30, 60, 90, 40, 80, 120)
c = outer(k,e)
@test c == ThirdOrderTensorValue{4,1,3}(10, 20, 30, 40, 20, 40, 60, 80, 30, 60, 90, 120)
@test_throws ArgumentError tr(c)

e = TensorValue{1,3}(10,20,30)
k = TensorValue{4,1}(1,2,3,4)
c = outer(e,k)
@test c == HighOrderTensorValue{Tuple{1,3,4,1}}(10, 20, 30, 20, 40, 60, 30, 60, 90, 40, 80, 120)
c = outer(k,e)
@test c == HighOrderTensorValue{Tuple{4,1,1,3}}(10, 20, 30, 40, 20, 40, 60, 80, 30, 60, 90, 120)
@test_throws ArgumentError tr(c)

e = VectorValue(10,20)
k = TensorValue(1,2,3,4)
@test tr(outer(e,k)) == VectorValue(50,110)

# Cross product

a = VectorValue(1,2,3)
b = VectorValue(3,1,6)
c = VectorValue(9,3,-5)
@test a × b == c

a = VectorValue(2,3)
b = VectorValue(5,2)
@test a × b == -11

a = VectorValue(1.0,5.0,-4.0)
b = VectorValue(6.0,2.0,3.0)
c = VectorValue(23.0,-27.0,-28.0)
@test cross(a, b) == c

a = VectorValue(4.0,1.0)
b = VectorValue(3.0,-2.0)
@test cross(a, b) == -11.0

a = VectorValue(4.0,1.0)
b = VectorValue(3.0,-2.0,1.0)
@test_throws ErrorException cross(a, b)

# Linear Algebra

t = TensorValue(10,2,30,4,5,6,70,8,9)

c = det(t)
@test c ≈ -8802.0
@test det(t) == det(TensorValue(get_array(t)))
@test inv(t) == inv(TensorValue(get_array(t)))

c = inv(t)
@test isa(c,TensorValue{3})

st = SymTensorValue(1:3...)
@test det(st) == det(TensorValue(get_array(st)))
@test inv(st) ≈  inv(TensorValue(get_array(st)))

st = SymTensorValue(9,8,7,5,4,1)
@test det(st) == det(TensorValue(get_array(st)))
@test inv(st) ≈  inv(TensorValue(get_array(st)))

qk = SkewSymTensorValue{1,Int}()
@test det(qk) == det(TensorValue(get_array(qk)))
@test inv(qk) ≈  inv(TensorValue(get_array(qk)))

qt = SymTracelessTensorValue(9,8,7,5,4)
@test det(qt) == det(TensorValue(get_array(qt)))
@test inv(qt) ≈  inv(TensorValue(get_array(qt)))

sk = SkewSymTensorValue{1,Int}()
@test det(sk) == det(TensorValue(get_array(sk)))
@test inv(sk) ≈  inv(TensorValue(get_array(sk)))

sk = SkewSymTensorValue(2.)
@test det(sk) == det(TensorValue(get_array(sk)))
@test inv(sk) ≈  inv(TensorValue(get_array(sk)))

sk = SkewSymTensorValue(1. : 3. ...)
@test det(sk) == det(TensorValue(get_array(sk)))
@test inv(sk) ==  SkewSymTensorValue(Inf,Inf,Inf)

sk = SkewSymTensorValue(1:6...)
@test det(sk) == det(TensorValue(get_array(sk)))
@test inv(sk) ≈  inv(TensorValue(get_array(sk)))

t = TensorValue(10)
@test det(t) == 10
@test inv(t) == TensorValue(1/10)

t = TensorValue(1,4,-1,1)
@test det(t) == det(TensorValue(get_array(t)))
@test inv(t) == inv(TensorValue(get_array(t)))

t = TensorValue(1:16...)
t += one(t)
@test det(t) == det(TensorValue(get_array(t)))
@test inv(t) == inv(TensorValue(get_array(t)))

q = SymTracelessTensorValue(1,2)
@test det(q) == det(TensorValue(get_array(q)))
@test inv(q) == SymTracelessTensorValue(inv(get_array(q)))

sk = SkewSymTensorValue(1,2,3)
@test det(sk) == 0

t = TensorValue{2,3}(1:6...)
@test_throws ErrorException det(t)
@test_throws ErrorException inv(t)

# Test eigen function
t = TensorValue(1.0,2.0,3.0,4.0)
result = eigen(t)
expected = eigen(get_array(t))
@test result.values ≈ expected.values
@test result.vectors ≈ expected.vectors

# Test different tensor types with eigen
st = SymTensorValue(9.0,8.0,7.0,5.0,4.0,1.0)
result_st = eigen(st)
expected_st = eigen(get_array(st))
@test result_st.values ≈ expected_st.values
@test result_st.vectors ≈ expected_st.vectors

qt = SymTracelessTensorValue(9.0,8.0,7.0,5.0,4.0)
result_qt = eigen(qt)
expected_qt = eigen(get_array(qt))
@test result_qt.values ≈ expected_qt.values
@test result_qt.vectors ≈ expected_qt.vectors

# Test 1x1 case
t1 = TensorValue(5.0)
result1 = eigen(t1)
expected1 = eigen(get_array(t1))
@test result1.values ≈ expected1.values
@test result1.vectors ≈ expected1.vectors

# Test 3x3 case
t3 = TensorValue(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,10.0)
result3 = eigen(t3)
expected3 = eigen(get_array(t3))
@test result3.values ≈ expected3.values
@test result3.vectors ≈ expected3.vectors

# Test error handling for non-square matrices
t_nonsquare = TensorValue{2,3}(1:6...)
@test_throws ErrorException eigen(t_nonsquare)

# Test sqrt function
t = TensorValue(4.0,0.0,0.0,4.0)
@test sqrt(t) == TensorValue(2.0,0.0,0.0,2.0)
t = SymTensorValue{2}(1,0,0)
@test sqrt(t) == TensorValue(1.0,0.0,0.0,0.0)
t = SymTracelessTensorValue{2,ComplexF64}(1,0)
@test sqrt(t) == TensorValue(1.0, 0.0, 0.0, -im)
t = SymTensorValue{2}(-1,0,0)
@test_throws DomainError sqrt(t)

# Measure

a = VectorValue(1,2,3)
c = meas(a)
@test c ≈ 3.7416573867739413

t = TensorValue(10,2,30,4,5,6,70,8,9)
c = meas(t)
@test c ≈ 8802.0

st = SymTensorValue(1,2,3,5,6,9)
@test meas(st) == meas(TensorValue(get_array(st)))

qt = SymTracelessTensorValue(1,2,3,5,6)
@test meas(qt) == meas(TensorValue(get_array(qt)))

sk = SkewSymTensorValue(1,2,3)
@test meas(sk) == meas(TensorValue(get_array(sk)))

v = TensorValue{1,2}(10,20)
@test meas(v) == sqrt(500)

v = TensorValue{2,3}(1,0,0,1,0,0)
@test meas(v) ≈ 1.0

v = TensorValue{2,3}(1,0,0,1,1,0)
@test meas(v) ≈ sqrt(2)

# Broadcasted operations

a = VectorValue(1,2,3)
b = VectorValue(2,1,6)

A = Array{VectorValue{3,Int},2}(undef,(4,5))
A .= a
@test A == fill(a,(4,5))

C = .- A

r = VectorValue(-1,-2,-3)
R = fill(r,(4,5))
@test C == R

B = fill(b,(4,1))

C = A .+ B

r = VectorValue(3,3,9)
R = fill(r,(4,5))
@test C == R

C = A .- B

r = VectorValue(-1,1,-3)
R = fill(r,(4,5))
@test C == R

# conj

v = VectorValue(1,0)
@test v == v'
@test conj(v) == v'

# tr

t = TensorValue(1,2,3,4)
@test tr(t) == 5

t = TensorValue(1,2,3,4,5,6,7,8,9)
@test tr(t) == 15

st = SymTensorValue(1,2,3,5,6,9)
@test tr(st) == tr(TensorValue(get_array(st)))

qt = SymTracelessTensorValue(1,2,3,5,6)
@test tr(qt) == tr(TensorValue(get_array(qt)))

sk = SkewSymTensorValue(1,2,3)
@test tr(sk) == tr(TensorValue(get_array(sk)))

t3 = ThirdOrderTensorValue(1:8...)
@test tr(t3) == tr(MultiValue(get_array(t3 )))

t23 = TensorValue{2,3}(1:6...)
@test_throws ArgumentError tr(t23)

t4 = HighOrderTensorValue(Val(4), 1:16...)
@test_throws ArgumentError  tr(t4)

# symmetric and skew-symmetric parts

@test get_array(symmetric_part(t)) == get_array(TensorValue(1.0, 3.0, 5.0, 3.0, 5.0, 7.0, 5.0, 7.0, 9.0))
@test symmetric_part(st) == symmetric_part(TensorValue(get_array(st)))
@test symmetric_part(st) === st
@test symmetric_part(qt) == symmetric_part(TensorValue(get_array(qt)))
@test symmetric_part(qt) === qt
ss = symmetric_part(sk)
@test zero(ss) == ss == symmetric_part(TensorValue(get_array(sk)))

@test skew_symmetric_part(sk) == SkewSymTensorValue(get_array(sk))
@test sk === skew_symmetric_part(sk)
skt = .5(t - t')
@test skew_symmetric_part(skt) == SkewSymTensorValue(get_array(skt))
skt = .5(st - st')
sk = skew_symmetric_part(skt)
@test zero(sk) == sk == SkewSymTensorValue(get_array(skt))
skt = .5(qt - qt')
sk = skew_symmetric_part(skt)
@test zero(sk) == sk == SkewSymTensorValue(get_array(skt))

# adjoint and transpose

a = TensorValue(1,2,3,4)
b = a'
@test adjoint(a) == b
@test b == TensorValue(1,3,2,4)
@test a⋅b == TensorValue(10,14,14,20)

a = TensorValue(1+0im,2-im,3,4+2im)
b = a'
@test adjoint(a) == b
@test b == TensorValue(1,3,2+im,4-2im)
@test a⋅b == TensorValue(10,14+5im,14-5im,25)

a = TensorValue(1,2,3,4)
b = a'
@test transpose(a) == b
@test transpose(a) == adjoint(a)
@test b == TensorValue(1,3,2,4)
@test a⋅b == TensorValue(10,14,14,20)

sa = SymTensorValue(1,2,3,5,6,9)
sb = sa'
@test adjoint(sa) == sb
@test sb == SymTensorValue(1,2,3,5,6,9)
@test sa⋅sb == TensorValue(get_array(sa))⋅TensorValue(get_array(sb))

sa = SymTensorValue(1+im,2,3,5,6,9)
sb = SymTensorValue(1-im,2,3,5,6,9)
@test sa' == sb
@test adjoint(sa) == sb
@test transpose(sa) === sa
@test sa⋅sb == TensorValue(get_array(sa))⋅TensorValue(get_array(sb))

sa = SymTracelessTensorValue(1+im,2,3,5,6)
sb = SymTracelessTensorValue(1-im,2,3,5,6)
@test sa' == sb
@test adjoint(sa) == sb
@test transpose(sa) === sa
@test sa⋅sb == TensorValue(get_array(sa))⋅TensorValue(get_array(sb))

sa = SkewSymTensorValue(1+im,2,3)
sb = SkewSymTensorValue(-1+im,-2,-3)
@test sa' == sb
@test adjoint(sa) == sb
@test transpose(sa) === -sa
@test sa⋅sb == TensorValue(get_array(sa))⋅TensorValue(get_array(sb))

u = VectorValue()
v = VectorValue{0,ComplexF64}()
@test norm(u) == 0
@test norm(v) == 0. + 0im

u = VectorValue(1.0,2.0)
v = VectorValue(2.0,3.0)
@test dot(u,v) ≈ inner(u,v)
@test norm(u) ≈ sqrt(inner(u,u))

u = VectorValue(1.0,im)
@test norm(u) ≈ sqrt(2)

a = TensorValue(1,2,3,4)
a2 = TensorValue(1,im,1,im)
@test norm(a) ≈ sqrt(inner(a,a))
@test norm(a2) ≈ 2

u = VectorValue(3.0,4.0)
@test normalize(u) ≈ VectorValue(0.6,0.8)

a = VectorValue(1.0,2.0)
b = VectorValue(2.0,3.0)
@test [a,b] ≈ [a,b]

a = VectorValue{0,Int}()
@test a ≈ a

λ = 1
μ = 1
ε = SymTensorValue(1,2,3)
σ = λ*tr(ε)*one(ε) + 2*μ*ε
@test isa(σ,SymTensorValue)
@test (σ ⊙ ε) == 52
#@test σ:ε == 52

I = one(SymFourthOrderTensorValue{2,Int})
@test I[1,1,1,1] == 1
@test I[1,2,1,2] == 0.5
@test I[2,1,1,2] == 0.5
@test I[2,2,2,2] == 1

@test I ⊙ ε == ε
@test ε ⊙ I == ε
@test I ⋅² ε == ε
@test ε ⋅² I == ε

a = TensorValue(1,2,3,4)
b = I ⊙ a
c = I ⋅² a
@test c == b == symmetric_part(a)


σ1 = λ*tr(ε)*one(ε) + 2*μ*ε
C = 2*μ*one(ε⊗ε) + λ*one(ε)⊗one(ε)
σ2 = C ⊙ ε
σ3 = C ⋅² ε
@test σ1 == σ2 == σ3

# third order contraction testing

# >1000 prime numbers
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, 8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, 8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, 8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831, 8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, 9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, 9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, 9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733, 9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973]

a = reshape(primes[4:30],(3,3,3))
a_tensor = ThirdOrderTensorValue(a...)
left_contraction = Matrix(get_array(VectorValue(1,2,-1) ⋅ a_tensor))
left_contraction_array = a[1,:,:] + 2*a[2,:,:] + -1*a[3,:,:]
@test left_contraction == left_contraction_array
right_contraction = Matrix(get_array(a_tensor ⋅ VectorValue(1,2,-1)))
right_contraction_array = a[:,:,1] + 2*a[:,:,2] + -1*a[:,:,3]
@test right_contraction == right_contraction_array

a = ThirdOrderTensorValue{3,3,3}(primes[1:27]...)
b = TensorValue(primes[28:36]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)
c = MultiValue(SArray(b)) ⋅² MultiValue(SArray(a))
@test SArray(c) == SArray(b ⋅² a)

a = ThirdOrderTensorValue{3,3,3}(primes[1:27]...)
b = SymTensorValue(primes[28:33]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)
c = MultiValue(SArray(b)) ⋅² MultiValue(SArray(a))
@test SArray(c) == SArray(b ⋅² a)

a = ThirdOrderTensorValue{3,3,3}(primes[1:27]...)
b = SymTracelessTensorValue(primes[28:32]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)
c = MultiValue(SArray(b)) ⋅² MultiValue(SArray(a))
@test SArray(c) == SArray(b ⋅² a)

a = ThirdOrderTensorValue{3,3,3}(primes[1:27]...)
b = SkewSymTensorValue(primes[28:30]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)
c = MultiValue(SArray(b)) ⋅² MultiValue(SArray(a))
@test SArray(c) == SArray(b ⋅² a)


# double Contractions w/ products

v  = VectorValue(primes[1:2]...)
t1 = TensorValue(primes[3:6]...)
t2 = TensorValue(primes[7:15]...)
s4 = SymFourthOrderTensorValue(primes[7:15]...)
@test_throws DimensionMismatch double_contraction(t1,v)
@test_throws DimensionMismatch double_contraction(t1,t2)

Sym4TensorIndexing = [1111, 1121, 1131, 1122, 1132, 1133, 2111, 2121, 2131, 2122, 2132, 2133,
                      3111, 3121, 3131, 3122, 3132, 3133, 2211, 2221, 2231, 2222, 2232, 2233,
                      2311, 2321, 2331, 2322, 2332, 2333, 3311, 3321, 3331, 3322, 3332, 3333]
test1 = SymFourthOrderTensorValue(primes[1:36]...)
test2 = SymFourthOrderTensorValue(primes[37:72]...)
result = Int[]
for off_index in Sym4TensorIndexing
  i = parse(Int,string(off_index)[1]); j = parse(Int,string(off_index)[2]);
  m = parse(Int,string(off_index)[3]); p = parse(Int,string(off_index)[4]);
  inner = 0
  for k in [1,2,3]
    for l in [1,2,3]
      inner += test1[i,j,k,l]*test2[k,l,m,p]
    end
  end
  append!(result,inner)
end
@test result == [(SymFourthOrderTensorValue(primes[1:36]...) ⋅² SymFourthOrderTensorValue(primes[37:72]...))...]

a = SymFourthOrderTensorValue{2}(primes[1:9]...)
b = SymFourthOrderTensorValue{2}(primes[10:18]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

a = SymFourthOrderTensorValue{3}(primes[1:36]...)
b = SymFourthOrderTensorValue{3}(primes[37:72]...)
#c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
c = HighOrderTensorValue{Tuple{3, 3, 3, 3}}([16058 57615 111407; 57615 167523 227169; 111407 227169 302217;;; 16442 58901 113877; 58901 171169 232067; 113877 232067 308771;;; 16994 60717 117341; 60717 176277 238911; 117341 238911 317939;;;; 16442 58901 113877; 58901 171169 232067; 113877 232067 308771;;; 17206 61551 119023; 61551 178839 242445; 119023 242445 322641;;; 17448 62647 121287; 62647 182399 247393; 121287 247393 329129;;;; 16994 60717 117341; 60717 176277 238911; 117341 238911 317939;;; 17448 62647 121287; 62647 182399 247393; 121287 247393 329129;;; 17818 64183 124391; 64183 187211 254049; 124391 254049 337921])
@test SArray(c) == SArray(a ⋅² b)

a = SymFourthOrderTensorValue{4}(primes[1:100]...)
b = SymFourthOrderTensorValue{4}(primes[101:200]...)
#c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
c = HighOrderTensorValue{Tuple{4, 4, 4, 4}, Int64}([188413 712136 1304296 2026734; 712136 2748622 3488732 4312968; 1304296 3488732 5129196 5935236; 2026734 4312968 5935236 6749138;;; 189575 717188 1313908 2041950; 717188 2769358 3515192 4345748; 1313908 3515192 5168440 5980840; 2041950 4345748 5980840 6800930;;; 190827 722090 1322958 2056060; 722090 2788564 3539610 4375934; 1322958 3539610 5204362 6022438; 2056060 4375934 6022438 6848208;;; 191791 725884 1330012 2067110; 725884 2803590 3558720 4399596; 1330012 3558720 5232580 6055128; 2067110 4399596 6055128 6885390;;;; 189575 717188 1313908 2041950; 717188 2769358 3515192 4345748; 1313908 3515192 5168440 5980840; 2041950 4345748 5980840 6800930;;; 193313 732352 1342192 2086054; 732352 2829214 3591684 4440152; 1342192 3591684 5281152 6111468; 2086054 4440152 6111468 6949414;;; 194655 738228 1353364 2103654; 738228 2853238 3622320 4478092; 1353364 3622320 5326516 6164152; 2103654 4478092 6164152 7009246;;; 196095 743898 1363934 2120148; 743898 2875580 3650802 4513262; 1363934 3650802 5368478 6212798; 2120148 4513262 6212798 7064468;;;; 190827 722090 1322958 2056060; 722090 2788564 3539610 4375934; 1322958 3539610 5204362 6022438; 2056060 4375934 6022438 6848208;;; 194655 738228 1353364 2103654; 738228 2853238 3622320 4478092; 1353364 3622320 5326516 6164152; 2103654 4478092 6164152 7009246;;; 197749 750626 1376534 2139800; 750626 2902328 3684882 4555438; 1376534 3684882 5418786 6271086; 2139800 4555438 6271086 7130732;;; 199085 756340 1387284 2156610; 756340 2925218 3714192 4591660; 1387284 3714192 5462008 6321248; 2156610 4591660 6321248 7187638;;;; 191791 725884 1330012 2067110; 725884 2803590 3558720 4399596; 1330012 3558720 5232580 6055128; 2067110 4399596 6055128 6885390;;; 196095 743898 1363934 2120148; 743898 2875580 3650802 4513262; 1363934 3650802 5368478 6212798; 2120148 4513262 6212798 7064468;;; 199085 756340 1387284 2156610; 756340 2925218 3714192 4591660; 1387284 3714192 5462008 6321248; 2156610 4591660 6321248 7187638;;; 200937 763194 1399646 2175840; 763194 2951440 3747226 4632710; 1399646 3747226 5510550 6377398; 2175840 4632710 6377398 7251424])
@test SArray(c) == SArray(a ⋅² b)

a = ThirdOrderTensorValue{3,3,3}(primes[1:27]...)
b = SymFourthOrderTensorValue{3}(primes[28:63]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_ilm = b_ij*c_jlm
a = TensorValue{3,3}(primes[1:9]...)
b = ThirdOrderTensorValue{3,3,3}(primes[10:36]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_il = b_ijk*c_jkl
a = ThirdOrderTensorValue{3,3,3}(primes[1:27]...)
b = TensorValue{3,3}(primes[28:36]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_il = b_ijk*c_jkl
a = ThirdOrderTensorValue{3,2,2}(primes[1:12]...)
b = ThirdOrderTensorValue{2,2,1}(primes[13:16]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_kl = b_ij*c_ijkl
a = SymTensorValue{3}(primes[1:6]...)
b = SymFourthOrderTensorValue(primes[7:42]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_ij = b_ijkl*c_kl
a = SymFourthOrderTensorValue(primes[1:36]...)
b = SymTensorValue{3}(primes[37:42]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_ij = b_ijkl*c_kl
c = MultiValue(SArray(b)) ⋅² MultiValue(SArray(a))
@test SArray(c) == SArray(b ⋅² a)

# a_k = b_ij*c_ijk
a = TensorValue{3,2}(primes[1:6]...)
b = ThirdOrderTensorValue{3,2,4}(primes[7:30]...)
c = MultiValue(SArray(a)) ⋅² MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅² b)

# a_il = b_ij*c_jl
a = TensorValue(primes[1:9]...)
b = TensorValue(primes[10:18]...)
c = MultiValue(SArray(a)) ⋅ MultiValue(SArray(b))
@test SArray(c) == SArray(a ⋅ b)
c = MultiValue(SArray(b)) ⋅ MultiValue(SArray(a))
@test SArray(c) == SArray(b ⋅ a)

# product of congruence

t = TensorValue{2,3}(1:6...)
t2 = TensorValue{2,2}(1:4...)
s = SymTensorValue{2}(1,2,3)
st = SymTracelessTensorValue{2}(1,2)
sk = SkewSymTensorValue{2}(1)

r = congruent_prod(t2,t)
@test r isa TensorValue{3,3}
@test r == t'⋅t2⋅t

r = congruent_prod(s,t)
@test r isa SymTensorValue{3}
@test get_array(r) == get_array(t'⋅s⋅t)

r = congruent_prod(sk,t)
@test r isa SkewSymTensorValue{3}
@test get_array(r) == get_array(t'⋅sk⋅t)

r = congruent_prod(st,t)
@test r isa SymTensorValue{3}
@test get_array(r) == get_array(t'⋅st⋅t)

@test_throws ErrorException congruent_prod(t,t2) # incompatible size

# Complex

a = 1.0 + 3.0*im
b = 4.0 - 3.0*im
@test outer(a,b) == a*b
@test inner(a,b) == a*b

v1 = VectorValue(1+1im)
v2 = VectorValue(1)
@test real(v1) == v2 && eltype(real(v1)) == eltype(v2)
@test imag(v1) == v2 && eltype(imag(v1)) == eltype(v2)
v2 = VectorValue(1-1im)
@test conj(v1) == v2 && eltype(conj(v1)) == eltype(v2)

v1 = VectorValue(1+1im, 1+1im)
v2 = VectorValue(1, 1)
@test real(v1) == v2 && eltype(real(v1)) == eltype(v2)
@test imag(v1) == v2 && eltype(imag(v1)) == eltype(v2)
v2 = VectorValue(1-1im, 1-1im)
@test conj(v1) == v2 && eltype(conj(v1)) == eltype(v2)

v1 = VectorValue(1+1im, 1+1im, 1+1im)
v2 = VectorValue(1, 1, 1)
@test real(v1) == v2 && eltype(real(v1)) == eltype(v2)
@test imag(v1) == v2 && eltype(imag(v1)) == eltype(v2)
v2 = VectorValue(1-1im, 1-1im, 1-1im)
@test conj(v1) == v2 && eltype(conj(v1)) == eltype(v2)

v1 = TensorValue(1+1im, 1+1im, 1+1im, 1+1im)
v2 = TensorValue(1, 1, 1, 1)
@test real(v1) == v2 && eltype(real(v1)) == eltype(v2)
@test imag(v1) == v2 && eltype(imag(v1)) == eltype(v2)
v2 = TensorValue(1-1im, 1-1im, 1-1im, 1-1im)
@test conj(v1) == v2 && eltype(conj(v1)) == eltype(v2)

v1 = SymTracelessTensorValue(1+1im, 1+1im)
v2 = SymTracelessTensorValue(1, 1)
@test real(v1) == v2 && eltype(real(v1)) == eltype(v2)
@test imag(v1) == v2 && eltype(imag(v1)) == eltype(v2)
v2 = SymTracelessTensorValue(1-1im, 1-1im)
@test conj(v1) == v2 && eltype(conj(v1)) == eltype(v2)

v1 = SkewSymTensorValue(1+1im)
v2 = SkewSymTensorValue(1)
@test real(v1) == v2 && eltype(real(v1)) == eltype(v2)
@test imag(v1) == v2 && eltype(imag(v1)) == eltype(v2)
v2 = SkewSymTensorValue(1-1im)
@test conj(v1) == v2 && eltype(conj(v1)) == eltype(v2)

#_eltype
@test _eltype(real,Tuple{},VectorValue(1.0,2.0))==Union{}
@test _eltype(real,Tuple{},VectorValue(1.0,2.0),VectorValue(2.0,3.0),VectorValue(3.0,4.0))==Union{}
a=VectorValue(4.0,5.0)
@test _eltype(real,a,VectorValue(1.0,2.0+im),VectorValue(2.0,3.0),VectorValue(3.0,4.0))==eltype(a)

# Broadcast
a = VectorValue(1,2,3)
b = VectorValue(1.,2.,3.)
c = a .* b
@test isa(c,VectorValue)
@test c.data == map(*,a.data,b.data)

a = TensorValue(1,2,3,4)
b = TensorValue(1.,2.,3.,4.)
c = a .* b
@test isa(c,TensorValue)
@test c.data == map(*,a.data,b.data)

@test diag(a) == VectorValue(1,4)


a = SymTensorValue(1,2,4)
b = SymTensorValue(1.,2.,4.)
c = a .* b
@test isa(c,SymTensorValue)
@test c.data == map(*,a.data,b.data)

@test diag(a) == VectorValue(1,4)


# Componant wise operations on sym. traceless tensors yield sym. tensors
a = SymTracelessTensorValue(1,2)
b = SymTracelessTensorValue(1.,2.)
c = a .* b
@test isa(c,SymTensorValue)
@test c.data == map(*,a.data,b.data)

a = SkewSymTensorValue(1,2,3)
@test diag(a) == zero(VectorValue(1:3...))

# Operation involving empty tensors (check type promotion)

v0  = VectorValue{0,Float64}()
v2  = VectorValue{2,Float64}(0,0)
t00 = TensorValue{0,0,Float32}()
t20 = TensorValue{2,0,Float32}()
t02 = TensorValue{0,2,Float32}()
t23 = TensorValue{2,3,Float32}(0,0,0,0,0,0)
st0 = SymTracelessTensorValue{0,ComplexF64}()
sk0 = SkewSymTensorValue{0,ComplexF64}()
t210= ThirdOrderTensorValue{2,1,0,ComplexF16}()
t012= ThirdOrderTensorValue{0,1,2,ComplexF16}()
t123= ThirdOrderTensorValue{1,2,3,Float64}(0,0,0,0,0,0)
f0  = SymFourthOrderTensorValue{0,ComplexF16}()

@test v0  ⋅ t02   === VectorValue{2,Float64}(0,0)
@test t20 ⋅ v0    === VectorValue{2,Float64}(0,0)
@test v0  ⋅ st0   === VectorValue{0,ComplexF64}()
@test st0 ⋅ v0    === VectorValue{0,ComplexF64}()
@test v0  ⋅ sk0   === VectorValue{0,ComplexF64}()
@test sk0 ⋅ v0    === VectorValue{0,ComplexF64}()
@test t20 ⋅ t02   === TensorValue{2,2,Float32}(0,0,0,0)
@test t02 ⋅ t20   === TensorValue{0,0,Float32}()
@test v0   ⋅ t012 === TensorValue{1,2,ComplexF64}(0,0)
@test t210 ⋅ v0   === TensorValue{2,1,ComplexF64}(0,0)
@test v2   ⋅ t210 === TensorValue{1,0,ComplexF64}()
@test t012 ⋅ v2   === TensorValue{0,1,ComplexF64}()
@test t012 ⋅ t23  === ThirdOrderTensorValue{0,1,3,ComplexF32}()
@test t02  ⋅ t210 === ThirdOrderTensorValue{0,1,0,ComplexF32}()

@test sk0 ⊙ sk0   === zero(ComplexF64)
@test st0 ⊙ st0   === zero(ComplexF64)
@test st0 ⊙ t00   === zero(ComplexF64)
@test t00 ⊙ st0   === zero(ComplexF64)
@test st0 ⊙ sk0   === zero(ComplexF64)
@test sk0 ⊙ st0   === zero(ComplexF64)
@test sk0 ⊙ t00   === zero(ComplexF64)
@test t00 ⊙ sk0   === zero(ComplexF64)
@test t00 ⊙ t00   === zero(Float32)
@test t00 ⊙ f0    === SymTensorValue{0,ComplexF32,0}()
@test f0  ⊙ t00   === SymTensorValue{0,ComplexF32,0}()
@test f0  ⊙ f0    === zero(ComplexF16)

t12  = TensorValue{1,2,Float32}(0,0)
t10  = TensorValue{1,0,Float32}()
st2  = SymTracelessTensorValue{2,Float32}(0,0)
sk2  = SkewSymTensorValue{2,Float32}(0)
t022 = ThirdOrderTensorValue{0,2,2,ComplexF16}()
t220 = ThirdOrderTensorValue{2,2,0,ComplexF16}()
t200 = ThirdOrderTensorValue{2,0,0,Float64}()
t002 = ThirdOrderTensorValue{0,0,2,Float64}()

@test t00  ⋅² t00  === zero(Float32)
@test t00  ⋅² f0   === SymTensorValue{0,ComplexF32,0}()
@test f0   ⋅² t00  === SymTensorValue{0,ComplexF32,0}()
@test st0  ⋅² f0   === SymTensorValue{0,ComplexF64,0}()
@test f0   ⋅² st0  === SymTensorValue{0,ComplexF64,0}()
@test sk0  ⋅² f0   === SymTensorValue{0,ComplexF64,0}()
@test f0   ⋅² sk0  === SymTensorValue{0,ComplexF64,0}()
@test t012 ⋅² t12  === VectorValue{0,ComplexF32}()
@test t210 ⋅² t10  === VectorValue{2,ComplexF32}(0,0)
@test t200 ⋅² st0  === VectorValue{2,ComplexF64}(0,0)
@test t022 ⋅² st2  === VectorValue{0,ComplexF32}()
@test t022 ⋅² sk2  === VectorValue{0,ComplexF32}()
@test st0  ⋅² t002 === VectorValue{2,ComplexF64}(0,0)
@test st2  ⋅² t220 === VectorValue{0,ComplexF32}()
@test sk2  ⋅² t220 === VectorValue{0,ComplexF32}()
@test t022 ⋅² t220 === TensorValue{0,0,ComplexF16}()
@test t200 ⋅² t002 === TensorValue{2,2,Float64}(0,0,0,0)
@test f0   ⋅² f0   === SymFourthOrderTensorValue{0,ComplexF16}()

@test 1im ⊗ v0  === VectorValue{0, ComplexF64}()
@test v0  ⊗ 1im === VectorValue{0, ComplexF64}()
@test v0  ⊗ v2  === TensorValue{0,2,Float64}()
@test v2  ⊗ v0  === TensorValue{2,0,Float64}()
@test v0  ⊗ t12 === ThirdOrderTensorValue{0,1,2,Float64}()
@test t00 ⊗ v2  === ThirdOrderTensorValue{0,0,2,Float64}()
@test st0 ⊗ st0 === SymFourthOrderTensorValue{0, ComplexF64}()

t10c = TensorValue{1,0,ComplexF32}()
t01c = TensorValue{0,1,ComplexF32}()

@test conj(st0) === st0
@test real(st0) === SymTracelessTensorValue{0,Float64}()
@test imag(st0) === SymTracelessTensorValue{0,Float64}()
@test conj(sk0) === sk0
@test real(sk0) === SkewSymTensorValue{0,Float64}()
@test imag(sk0) === SkewSymTensorValue{0,Float64}()
@test tr(t00) === zero(Float32)
@test tr(t220) === VectorValue{0,ComplexF16}()
@test adjoint(t00) === t00
@test adjoint(t01c) === t10c
@test symmetric_part(t00) === SymTensorValue{0,Float32}()

# permutedims

# permutedims

v = VectorValue(1,2,3)
@test permutedims(v, (1,)) == v
@test isa(permutedims(v, (1,)), VectorValue{3,Int})

A = TensorValue{2,3}(1,2,3,4,5,6)
@test permutedims(A, (1,2)) == A
@test isa(permutedims(A, (1,2)), TensorValue{2,3,Int})
At = permutedims(A, (2,1))
@test isa(At, TensorValue{3,2,Int})
@test At == transpose(A)

S = SymTensorValue{3}(1,2,3,4,5,6)
@test permutedims(S, (2,1)) == TensorValue{3,3,Int}(get_array(S)')

T3 = ThirdOrderTensorValue{2,3,4}(1:24...)
@test permutedims(T3, (1,2,3)) == HighOrderTensorValue{Tuple{2,3,4},Int,3}(Tuple(get_array(T3)))

P213 = permutedims(T3, (2,1,3))
@test isa(P213, HighOrderTensorValue{Tuple{3,2,4},Int,3})
@test all(P213[j,i,k] == T3[i,j,k] for i in 1:2, j in 1:3, k in 1:4)

P312 = permutedims(T3, (3,1,2))
@test isa(P312, HighOrderTensorValue{Tuple{4,2,3},Int,3})
@test all(P312[k,i,j] == T3[i,j,k] for i in 1:2, j in 1:3, k in 1:4)

r1 = tensor_contraction(T3, VectorValue(1,2,3), (2,), (1,))
r2 = tensor_contraction(P213, VectorValue(1,2,3), (1,), (1,))
@test r1 == r2

@test_throws ArgumentError permutedims(A, (1,3))
@test_throws ArgumentError permutedims(A, (1,1))

# tensor_contraction (self-contraction)

# trace as special case
A2 = TensorValue{3,3}(1:9...)
@test tensor_contraction(A2, (1,), (2,)) == tr(A2)
@test isa(tensor_contraction(A2, (1,), (2,)), Int)

# non-square matrix: contract unique index pair → scalar
A23 = TensorValue{2,2}(1,2,3,4)
@test tensor_contraction(A23, (1,), (2,)) == tr(A23)

# partial trace of a 4th-order tensor: C[j,l] = Σ_i A[i,j,i,l]
A4 = HighOrderTensorValue{Tuple{3,2,3,4}}(reshape(1:72, 3,2,3,4))
C  = tensor_contraction(A4, (1,), (3,))
@test isa(C, TensorValue{2,4,Int})
@test all(C[j,l] == sum(A4[i,j,i,l] for i in 1:3) for j in 1:2, l in 1:4)

# double self-contraction: scalar
A4b = HighOrderTensorValue{Tuple{2,3,2,3}}(reshape(1:36, 2,3,2,3))
s   = tensor_contraction(A4b, (1,2), (3,4))
@test isa(s, Int)
@test s == sum(A4b[i,j,i,j] for i in 1:2, j in 1:3)

# consistency: tr(a) == tensor_contraction(a, (1,), (2,))
st = SymTensorValue{3}(1,2,3,4,5,6)
@test tensor_contraction(st, (1,), (2,)) == tr(st)

# integer signatures
@test tensor_contraction(A2, 1, 2) == tr(A2)
@test tensor_contraction(A2, VectorValue(1,2,3), 2, 1) == tensor_contraction(A2, VectorValue(1,2,3), (2,), (1,))

# error cases
@test_throws ArgumentError tensor_contraction(A2, (1,), (1,))          # i∩j ≠ ∅
@test_throws ArgumentError tensor_contraction(A2, (3,), (1,))          # out of range
@test_throws DimensionMismatch tensor_contraction(TensorValue{2,3}(1:6...), (1,), (2,))

# tensor_contraction (two tensors)

A = TensorValue{2,3}(1,2,3,4,5,6)
x = VectorValue(1,2,3)
r = tensor_contraction(A, x, (2,), (1,))
@test isa(r, VectorValue{2,Int})
@test r == A⋅x

y = VectorValue(1,2)
r2 = tensor_contraction(A, y, (1,), (1,))
@test isa(r2, VectorValue{3,Int})
ref = VectorValue(A[1,1]*y[1]+A[2,1]*y[2], A[1,2]*y[1]+A[2,2]*y[2], A[1,3]*y[1]+A[2,3]*y[2])
@test r2 == ref

u = VectorValue(1,2,3)
v = VectorValue(4,5,6)
@test tensor_contraction(u, v, (1,), (1,)) == u⋅v

M = TensorValue{2,3}(1,2,3,4,5,6)
N = TensorValue{2,3}(6,5,4,3,2,1)
@test tensor_contraction(M, N, (1,2), (1,2)) == contracted_product(Val(2), M, N)

t1 = TensorValue{2,3}(1,2,3,4,5,6)
t2 = TensorValue{2,3}(2,3,4,5,6,7)
@test tensor_contraction(t1, t2, (1,2), (1,2)) == inner(t1,t2)

M2 = TensorValue{2,3}(1,2,3,4,5,6)
N2 = TensorValue{3,2}(6,5,4,3,2,1)
@test tensor_contraction(M2, N2, (1,2), (2,1)) == sum(M2[i,j]*N2[j,i] for i in 1:2, j in 1:3)

a2 = VectorValue(1,2)
b3 = VectorValue(3,4,5)
@test tensor_contraction(a2, b3, (), ()) == outer(a2, b3)

T3 = ThirdOrderTensorValue{2,3,4}(1:24...)
w  = VectorValue(1,2,3)
r3 = tensor_contraction(T3, w, (2,), (1,))
@test isa(r3, TensorValue{2,4,Int})
ref3 = [sum(T3[i,k,j]*w[k] for k in 1:3) for i in 1:2, j in 1:4]
@test get_array(r3) == ref3

T3b = ThirdOrderTensorValue{3,2,3}(1:18...)
B   = TensorValue{3,3}(1:9...)
r4  = tensor_contraction(T3b, B, (1,3), (1,2))
@test isa(r4, VectorValue{2,Int})
ref4 = [sum(T3b[i,j,k]*B[i,k] for i in 1:3, k in 1:3) for j in 1:2]
@test r4 == VectorValue(ref4...)

@test_throws DimensionMismatch tensor_contraction(VectorValue(1,2), VectorValue(1,2,3), (1,), (1,))
@test_throws ArgumentError tensor_contraction(VectorValue(1,2), VectorValue(1,2), (3,), (1,))
@test_throws ArgumentError tensor_contraction(TensorValue{2,2}(1:4...), VectorValue(1,2), (1,1), (1,1))

end # module OperationsTests
