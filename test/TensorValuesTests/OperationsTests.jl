module OperationsTests

using Test: Error
using Test
using Gridap.TensorValues
using Gridap.Arrays
using LinearAlgebra
using Gridap.TensorValues: _eltype

# Comparison

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

a = VectorValue(1,2,3)
b = VectorValue(2,1,6)

@test a==a
@test a ≈ a
@test a!=b
@test [a,a] == [a,a]
@test [a,a] ≈ [a,a]

c = TensorValue(1,2,3,4)

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

v = VectorValue(1,2)
t = TensorValue(1,2,3,4)
s = SymTensorValue(1,2,3)
q = SymTracelessTensorValue(1,2)
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

# Matrix Division

a = VectorValue(1,2,3)

t = one(TensorValue{3,3,Int})
c = t\a
@test c == a

st = one(SymTensorValue{3,Int})
c = st\a
@test c == a

# Operations by a scalar

t = TensorValue(1,2,3,4,5,6,7,8,9)
st = SymTensorValue(1,2,3,5,6,9)
qt = SymTracelessTensorValue(1,2,3,5,6)
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

c = st ⋅ qt2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 78, 108, 18, 39, 54, -27, -52, -81)
@test c == r

c = qt2 ⋅ st
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 18, -27, 78, 39, -52, 108, 54, -81)
@test c == r

c = a ⋅ st
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,42)
@test c == r

c = a ⋅ qt
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,-3)
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
c = st ⊙ st2
@test isa(c,Int)
@test c == inner(TensorValue(get_array(st)),TensorValue(get_array(st2)))

c = inner(qt,qt2)
c = qt ⊙ qt2
@test isa(c,Int)
@test c == inner(TensorValue(get_array(qt)),TensorValue(get_array(qt2)))

t1 = TensorValue{2,3}(1:6...)
t2 = TensorValue{2,3}(10:15...)
@test inner(t1,t1) == 91
@test double_contraction(t1,t2) == inner(t1,t2)

@test_throws ErrorException inner(a,t)
@test_throws ErrorException inner(s,a)
@test_throws ErrorException inner(a,q)

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

st = SymTensorValue(9,8,7,5,4,1)
@test det(st) == det(TensorValue(get_array(st)))
@test inv(st) ≈  inv(TensorValue(get_array(st)))

qt = SymTracelessTensorValue(9,8,7,5,4)
@test det(qt) == det(TensorValue(get_array(qt)))
@test inv(qt) ≈  inv(TensorValue(get_array(qt)))

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

t = TensorValue{2,3}(1:6...)
@test_throws ErrorException det(t)
@test_throws ErrorException inv(t)

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

@test get_array(symmetric_part(t)) == get_array(TensorValue(1.0, 3.0, 5.0, 3.0, 5.0, 7.0, 5.0, 7.0, 9.0))
@test symmetric_part(st) == symmetric_part(TensorValue(get_array(st)))
@test symmetric_part(st) === st
@test symmetric_part(qt) == symmetric_part(TensorValue(get_array(qt)))
@test symmetric_part(qt) === qt

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

sa = SymTensorValue(1,2,3,5,6,9)
sb = sa'
@test transpose(sa) == sb
@test sb == SymTensorValue(1,2,3,5,6,9)
@test sa⋅sb == TensorValue(get_array(sa))⋅TensorValue(get_array(sb))

sa = SymTracelessTensorValue(1,2,3,5,6)
sb = sa'
@test adjoint(sa) == sb
@test sb == SymTracelessTensorValue(1,2,3,5,6)
@test sa⋅sb == TensorValue(get_array(sa))⋅TensorValue(get_array(sb))

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
a = reshape(Vector(1:27),(3,3,3))
a_tensor = ThirdOrderTensorValue(a...)
left_contraction = Matrix(get_array(VectorValue(1,2,-1) ⋅ a_tensor))
left_contraction_array = a[1,:,:] + 2*a[2,:,:] + -1*a[3,:,:]
@test left_contraction == left_contraction_array
right_contraction = Matrix(get_array(a_tensor ⋅ VectorValue(1,2,-1)))
right_contraction_array = a[:,:,1] + 2*a[:,:,2] + -1*a[:,:,3]
@test right_contraction == right_contraction_array

a = reshape(Vector(1:27),(3,3,3))
b = reshape(Vector(1:9),(3,3))
a_tensor = ThirdOrderTensorValue(a...)
b_tensor = TensorValue(b...)
odot_contraction = Vector(get_array(a_tensor ⋅² b_tensor))
odot_contraction_array = 1*a[:,1,1] + 4*a[:,1,2] + 7*a[:,1,3] + 2*a[:,2,1] +
  5*a[:,2,2] + 8*a[:,2,3] + 3*a[:,3,1] + 6*a[:,3,2] + 9*a[:,3,3]
@test odot_contraction == odot_contraction_array

a = reshape(Vector(1:27),(3,3,3))
a_tensor = ThirdOrderTensorValue(a...)
b_tensor = SymTensorValue((1:6)...)
b = Matrix(get_array(b_tensor))
odot_contraction = Vector(get_array(a_tensor ⋅² b_tensor))
odot_contraction_array = 1*a[:,1,1] + 2*a[:,1,2] + 3*a[:,1,3] + 2*a[:,2,1] +
  4*a[:,2,2] + 5*a[:,2,3] + 3*a[:,3,1] + 5*a[:,3,2] + 6*a[:,3,3]
@test odot_contraction == odot_contraction_array

a = reshape(Vector(1:27),(3,3,3))
a_tensor = ThirdOrderTensorValue(a...)
b_tensor = SymTracelessTensorValue((1:5)...)
b = Matrix(get_array(b_tensor))
odot_contraction = Vector(get_array(a_tensor ⋅² b_tensor))
odot_contraction_array = 1*a[:,1,1] + 2*a[:,1,2] + 3*a[:,1,3] + 2*a[:,2,1] +
  4*a[:,2,2] + 5*a[:,2,3] + 3*a[:,3,1] + 5*a[:,3,2] + (-5)*a[:,3,3]
@test odot_contraction == odot_contraction_array

# double Contractions w/ products

v  = VectorValue(1:2...)
t1 = TensorValue(1:4...)
t2 = TensorValue(1:9...)
s4 = SymFourthOrderTensorValue(1:9...)
@test_throws ErrorException double_contraction(t1,v)
@test_throws DimensionMismatch double_contraction(t1,t2)

Sym4TensorIndexing = [1111, 1121, 1131, 1122, 1132, 1133, 2111, 2121, 2131, 2122, 2132, 2133,
                      3111, 3121, 3131, 3122, 3132, 3133, 2211, 2221, 2231, 2222, 2232, 2233,
                      2311, 2321, 2331, 2322, 2332, 2333, 3311, 3321, 3331, 3322, 3332, 3333]
test1 = test2 = SymFourthOrderTensorValue(1:36...)
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
@test result == [(SymFourthOrderTensorValue(1:36...) ⋅² SymFourthOrderTensorValue(1:36...))...]

# generalised
test1 = test2 = SymFourthOrderTensorValue(1:9...)
res = zeros(2,2,2,2);
res[1,1,1,1]=1; res[1,1,1,2]=2; res[1,1,2,1]=2; res[1,1,2,2]=3;
res[1,2,1,1]=4; res[2,1,1,1]=4; res[1,2,1,2]=5; res[1,2,2,1]=5;
res[2,1,1,2]=5; res[2,1,2,1]=5; res[2,1,2,2]=6; res[1,2,2,2]=6;
res[2,2,2,2]=6; res[2,2,1,1]=7; res[2,2,2,1]=8; res[2,2,1,2]=8;
res[2,2,2,2]=9;
result = zeros(2,2,2,2)
for i in 1:2
  for j in 1:2
    for p in 1:2
      for m in 1:2
        result[i,j,p,m] = sum(res[i,j,k,l]*res[k,l,p,m] for k in 1:2 for l in 1:2)
      end
    end
  end
end
pass = true
for i in 1:2
  for j in 1:2
    for p in 1:2
      for m in 1:2
        if result[i,j,p,m] != (test1 ⋅² test2)[i,j,p,m]
          pass = false
        end
      end
    end
  end
end
@test pass

# a_ilm = b_ijk*c_jklm
vals = ones(3,3,3);
vals[1,:,:] .= [3 1 0; 1 4 0; 0 0 2];
vals[2,:,:] .= [3 0 2; 0 2 0; 2 0 1];
vals[3,:,:] .= [1 0 0; 0 2 1; 0 1 3];
t1 = ThirdOrderTensorValue(vals ...)
t2 = SymFourthOrderTensorValue(ones(36) ...)
for ci in CartesianIndices((3,3,3))
  @test (t1 ⋅² t2)[ci] == sum(t1[ci[1],j,k]*t2[j,k,ci[2],ci[3]] for j in 1:3 for k in 1:3)
end

# a_ilm = b_ij*c_jlm
vals = ones(3,3,3);
vals[1,:,:] .= [3 1 0
               1 4 0
               0 0 2];
vals[2,:,:] .= [3 0 2
                0 2 0
                2 0 1];
vals[3,:,:] .= [1 0 0
                0 2 1
                0 1 3];
t1 = TensorValue{3,3}(1:9...)
t2 = ThirdOrderTensorValue{3,3,3}(1:27...)
for ci in CartesianIndices((3,3,3))
  @test (t1 ⋅ t2)[ci] == sum(t1[ci[1],j]*t2[j,ci[2],ci[3]] for j in 1:3)
end

# a_il = b_ijk*c_jkl
vals = ones(3,3,3);
vals[1,:,:] .= [3 1 0
               1 4 0
               0 0 2];
vals[2,:,:] .= [3 0 2
                0 2 0
                2 0 1];
vals[3,:,:] .= [1 0 0
                0 2 1
                0 1 3];
t1 = ThirdOrderTensorValue(vals ...)
for ci in CartesianIndices((3,3))
  @test (t1 ⋅² t1)[ci] == sum(t1[ci[1],i,j]*t1[i,j,ci[2]] for i in 1:3 for j in 1:3)
end

# a_il = b_ijk*c_jkl
t1 = ThirdOrderTensorValue{3,2,2}(1:12...)
t2 = ThirdOrderTensorValue{2,2,1}(1:4...)
t1_double_t2 = t1 ⋅² t2
@test isa(t1_double_t2, TensorValue{3,1})
@test (t1 ⋅² t2)[1,1] == sum(t1[1,j,k] * t2[j,k,1] for j in 1:2 for k in 1:2)
@test (t1 ⋅² t2)[2,1] == sum(t1[2,j,k] * t2[j,k,1] for j in 1:2 for k in 1:2)
@test (t1 ⋅² t2)[3,1] == sum(t1[3,j,k] * t2[j,k,1] for j in 1:2 for k in 1:2)

# a_kl = b_ij*c_ijkl
t1 = SymTensorValue{3}(1:6...)
t2 = SymFourthOrderTensorValue(1:36 ...)
t1_double_t2 = t1 ⋅² t2
for ci in CartesianIndices((3,3))
  @test t1_double_t2[ci] == sum(t1[i,j]*t2[i,j,ci[1],ci[2]] for i in 1:3 for j in 1:3)
end

# a_ij = b_ijkl*c_kl
t1 = SymFourthOrderTensorValue(1:36...)
t2 = SymTensorValue{3}(1:6...)
t1_double_t2 = t1 ⋅² t2
for ci in CartesianIndices((3,3))
  @test t1_double_t2[ci] == sum(t1[ci[1],ci[2],i,j]*t2[i,j] for i in 1:3 for j in 1:3)
end

# a_ij = b_ijkl*c_kl
t2_double_t1 = t2 ⋅² t1
for ci in CartesianIndices((3,3))
  @test t2_double_t1[ci] == sum(t2[k,l]*t1[k,l,ci[1],ci[2]] for k in 1:3 for l in 1:3)
end

# a_k = b_ij*c_ijk
t1 = TensorValue{3,2}(1:6...)
t2 = ThirdOrderTensorValue{3,2,4}(1:24...)
t1_double_t2 = t1 ⋅² t2
for ci in CartesianIndices((4,))
  @test t1_double_t2[ci] == sum(t1[i,j]*t2[i,j,ci[1]] for i in 1:3 for j in 1:2)
end

# a_il = b_ij*c_jl
v1 = [1 2 3
      2 4 0
      3 0 5];
v2 = [1 1 2
      1 2 5
      2 5 3];
t1 = TensorValue(v1)
t2 = TensorValue(v2)
t1_dot_t2 = t1 ⋅ t2
for ci in CartesianIndices((3,3))
  @test t1_dot_t2[ci] == sum(t1[ci[1],i]*t2[i,ci[2]] for i in 1:3)
end
t2_dot_t1 = t2 ⋅ t1
for ci in CartesianIndices((3,3))
  @test t2_dot_t1[ci] == sum(t2[ci[1],i]*t1[i,ci[2]] for i in 1:3)
end

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

@test diag(a) == VectorValue(1,-1)

# Operation involving empty tensors (check type promotion)

v0  = VectorValue{0,Float64}()
v2  = VectorValue{2,Float64}(0,0)
t00 = TensorValue{0,0,Float32}()
t20 = TensorValue{2,0,Float32}()
t02 = TensorValue{0,2,Float32}()
t23 = TensorValue{2,3,Float32}(0,0,0,0,0,0)
st0 = SymTracelessTensorValue{0,ComplexF64}()
t210= ThirdOrderTensorValue{2,1,0,ComplexF16}()
t012= ThirdOrderTensorValue{0,1,2,ComplexF16}()
t123= ThirdOrderTensorValue{1,2,3,Float64}(0,0,0,0,0,0)
f0  = SymFourthOrderTensorValue{0,ComplexF16}()

@test v0  ⋅ t02   === VectorValue{2,Float64}(0,0)
@test t20 ⋅ v0    === VectorValue{2,Float64}(0,0)
@test v0  ⋅ st0   === VectorValue{0,ComplexF64}()
@test st0 ⋅ v0    === VectorValue{0,ComplexF64}()
@test t20 ⋅ t02   === TensorValue{2,2,Float32}(0,0,0,0)
@test t02 ⋅ t20   === TensorValue{0,0,Float32}()
@test v0   ⋅ t012 === TensorValue{1,2,ComplexF64}(0,0)
@test t210 ⋅ v0   === TensorValue{2,1,ComplexF64}(0,0)
@test v2   ⋅ t210 === TensorValue{1,0,ComplexF64}()
@test t012 ⋅ v2   === TensorValue{0,1,ComplexF64}()
@test t012 ⋅ t23  === ThirdOrderTensorValue{0,1,3,ComplexF32}()
@test t02  ⋅ t210 === ThirdOrderTensorValue{0,1,0,ComplexF32}()

@test st0 ⊙ st0   === zero(ComplexF64)
@test st0 ⊙ t00   === zero(ComplexF64)
@test t00 ⊙ st0   === zero(ComplexF64)
@test t00 ⊙ t00   === zero(Float32)
@test t00 ⊙ f0    === SymTensorValue{0,ComplexF32,0}()
@test f0  ⊙ t00   === SymTensorValue{0,ComplexF32,0}()
@test f0  ⊙ f0    === SymFourthOrderTensorValue{0, ComplexF16}()

t12  = TensorValue{1,2,Float32}(0,0)
t10  = TensorValue{1,0,Float32}()
st2  = SymTracelessTensorValue{2,Float32}(0,0)
t022 = ThirdOrderTensorValue{0,2,2,ComplexF16}()
t220 = ThirdOrderTensorValue{2,2,0,ComplexF16}()
t200 = ThirdOrderTensorValue{2,0,0,Float64}()
t002 = ThirdOrderTensorValue{0,0,2,Float64}()

@test t00  ⋅² t00  === zero(Float32)
@test t00  ⋅² f0   === SymTensorValue{0,ComplexF32,0}()
@test f0   ⋅² t00  === SymTensorValue{0,ComplexF32,0}()
@test st0  ⋅² f0   === SymTensorValue{0,ComplexF64,0}()
@test f0   ⋅² st0  === SymTensorValue{0,ComplexF64,0}()
@test t012 ⋅² t12  === VectorValue{0,ComplexF32}()
@test t210 ⋅² t10  === VectorValue{2,ComplexF32}(0,0)
@test t200 ⋅² st0  === VectorValue{2,ComplexF64}(0,0)
@test t022 ⋅² st2  === VectorValue{0,ComplexF32}()
@test st0  ⋅² t002 === VectorValue{2,ComplexF64}(0,0)
@test st2  ⋅² t220 === VectorValue{0,ComplexF32}()
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
@test tr(t00) === zero(Float32)
@test tr(t220) === VectorValue{0,ComplexF16}()
@test adjoint(t00) === t00
@test adjoint(t01c) === t10c
@test symmetric_part(t00) === SymTensorValue{0,Float32}()

end # module OperationsTests
