module OperationsTests

using Test
using Gridap.TensorValues
using Gridap.Arrays
using LinearAlgebra

# Comparison

a = VectorValue(1,2,3)
b = VectorValue(1,3,3)

@test (a < b) == true
@test (a <= b) == true
@test (a == b) == false
@test (a >= b) == false
@test (a > b) == false

a = VectorValue(1,2,3)
b = VectorValue(2,1,6)

@test a==a
@test a ≈ a
@test a!=b
@test [a,a] == [a,a]
@test [a,a] ≈ [a,a]

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

c = st + 2
@test isa(c,SymTensorValue{3})
r = SymTensorValue(3,4,5,7,8,11)
@test c == r

c = 2 * s4ot
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(2,0,0, 0,1,0, 0,0,2)
@test c == r

c = s4ot * 2
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(2,0,0, 0,1,0, 0,0,2)
@test c == r

c = c + 0
@test isa(c,SymFourthOrderTensorValue{2})
r = SymFourthOrderTensorValue(2,0,0, 0,1,0, 0,0,2)
@test c == r

# Dot product (simple contraction)

a = VectorValue(1,2,3)
b = VectorValue(2,1,6)

t = TensorValue(1,2,3,4,5,6,7,8,9)
s = TensorValue(9,8,3,4,5,6,7,2,1)
st = SymTensorValue(1,2,3,5,6,9)
st2 = SymTensorValue(9,6,5,3,2,1)

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

c = s ⋅ t
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(38,24,18,98,69,48,158,114,78)
@test c == r

c = st ⋅ st2
@test isa(c,TensorValue{3,3,Int})
r = TensorValue(36, 78, 108, 18, 39, 54, 12, 26, 36)
@test c == r

c = a ⋅ st
@test isa(c,VectorValue{3,Int})
r = VectorValue(14,30,42)
@test c == r

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

e = VectorValue(10,20)
k = TensorValue(1,2,3,4)
c = outer(e,k)
@test c == ThirdOrderTensorValue{2,2,2}(10, 20, 20, 40, 30, 60, 40, 80)

@test tr(c) == VectorValue(50,110)

# Cross product

a = VectorValue(1,2,3)
b = VectorValue(3,1,6)
c = VectorValue(9,3,-5)
@test a × b == c

a = VectorValue(2,3)
b = VectorValue(5,2)
@test a × b == -11

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
@test inv(st) == inv(TensorValue(get_array(st)))

t = TensorValue(10)
@test det(t) == 10
@test inv(t) == TensorValue(1/10)

t = TensorValue(1,4,-1,1)
@test det(t) == det(TensorValue(get_array(t)))
@test inv(t) == inv(TensorValue(get_array(t)))

# Measure

a = VectorValue(1,2,3)
c = meas(a)
@test c ≈ 3.7416573867739413

t = TensorValue(10,2,30,4,5,6,70,8,9)
c = meas(t)
@test c ≈ 8802.0

st = SymTensorValue(1,2,3,5,6,9)
@test meas(st) == meas(TensorValue(get_array(st)))

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

t = TensorValue(1,2,3,4)
@test tr(t) == 5

t = TensorValue(1,2,3,4,5,6,7,8,9)
@test tr(t) == 15

st = SymTensorValue(1,2,3,5,6,9)
@test tr(st) == tr(TensorValue(get_array(st)))

@test get_array(symmetric_part(t)) == get_array(TensorValue(1.0, 3.0, 5.0, 3.0, 5.0, 7.0, 5.0, 7.0, 9.0))
@test symmetric_part(st) == symmetric_part(TensorValue(get_array(st)))

a = TensorValue(1,2,3,4)
b = a'
@test adjoint(a) == b
@test b == TensorValue(1,3,2,4)
@test a⋅b == TensorValue(10,14,14,20)

a = TensorValue(1,2,3,4)
b = a'
@test transpose(a) == b
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

u = VectorValue(1.0,2.0)
v = VectorValue(2.0,3.0)
@test dot(u,v) ≈ inner(u,v)
@test norm(u) ≈ sqrt(inner(u,u))

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
#@test I : ε == ε

a = TensorValue(1,2,3,4)
b = I ⊙ a
@test b == symmetric_part(a)
#b = I : a
#@test b == symmetric_part(a)


σ1 = λ*tr(ε)*one(ε) + 2*μ*ε
C = 2*μ*one(ε⊗ε) + λ*one(ε)⊗one(ε)
σ2 = C ⊙ ε
@test σ1 == σ2
#σ2 = C : ε
#@test σ1 == σ2


end # module OperationsTests
