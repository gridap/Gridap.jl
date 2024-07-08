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

a = TensorValue(1,2,3,4)
b = SymTensorValue(5,6,7)

c = a + b
r = TensorValue(6,8,9,11)
@test c==r

a = SymTensorValue(5,6,7)
b = TensorValue(1,2,3,4)

c = a - b
r = TensorValue(4,4,3,3)
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

a = VectorValue(1.0,5.0,-4.0)
b = VectorValue(6.0,2.0,3.0)
c = VectorValue(23.0,-27.0,-28.0)
@test cross(a, b) == c

a = VectorValue(4.0,1.0)
b = VectorValue(3.0,-2.0)
@test cross(a, b) == -11.0

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

a = TensorValue(1,2,3,4)
@test norm(a) ≈ sqrt(inner(a,a))

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

# double Contractions w/ products
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
t2 = SymFourthOrderTensorValue(ones(36) ...)
v111 = sum(vals[1,j,k]*ones(3,3,3,3)[j,k,1,1] for j in 1:3 for k in 1:3);
v112 = sum(vals[1,j,k]*ones(3,3,3,3)[j,k,2,2] for j in 1:3 for k in 1:3);
v113 = sum(vals[1,j,k]*ones(3,3,3,3)[j,k,1,3] for j in 1:3 for k in 1:3);
v121 = sum(vals[1,j,k]*ones(3,3,3,3)[j,k,2,1] for j in 1:3 for k in 1:3);
v131 = sum(vals[1,j,k]*ones(3,3,3,3)[j,k,3,1] for j in 1:3 for k in 1:3);
v211 = sum(vals[2,j,k]*ones(3,3,3,3)[j,k,1,1] for j in 1:3 for k in 1:3);
v311 = sum(vals[3,j,k]*ones(3,3,3,3)[j,k,1,1] for j in 1:3 for k in 1:3);
@test v111 == (t1 ⋅² t2)[1,1,1]
@test v112 == (t1 ⋅² t2)[1,1,2]
@test v113 == (t1 ⋅² t2)[1,1,3]
@test v121 == (t1 ⋅² t2)[1,2,1]
@test v131 == (t1 ⋅² t2)[1,3,1]
@test v211 == (t1 ⋅² t2)[2,1,1]
@test v311 == (t1 ⋅² t2)[3,1,1]

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
t1 = TensorValue(ones(3,3)...)
t2 = ThirdOrderTensorValue(vals ...)
@test (t1 ⋅ t2)[1,1,1] == sum(ones(3,3)[1,j]*vals[j,1,1] for j in 1:3)
@test (t1 ⋅ t2)[3,2,1] == sum(ones(3,3)[3,j]*vals[j,2,1] for j in 1:3)
@test (t1 ⋅ t2)[1,1,2] == sum(ones(3,3)[1,j]*vals[j,1,2] for j in 1:3)
@test (t1 ⋅ t2)[1,1,3] == sum(ones(3,3)[1,j]*vals[j,1,3] for j in 1:3)
@test (t1 ⋅ t2)[1,2,2] == sum(ones(3,3)[1,j]*vals[j,2,2] for j in 1:3)
@test (t1 ⋅ t2)[1,2,3] == sum(ones(3,3)[1,j]*vals[j,2,3] for j in 1:3)
@test (t1 ⋅ t2)[1,3,3] == sum(ones(3,3)[1,j]*vals[j,3,3] for j in 1:3)

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
@test (t1 ⋅² t1)[1,1] == sum(vals[1,i,j] .* vals[i,j,1] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[2,1] == sum(vals[2,i,j] .* vals[i,j,1] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[3,1] == sum(vals[3,i,j] .* vals[i,j,1] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[1,2] == sum(vals[1,i,j] .* vals[i,j,2] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[2,2] == sum(vals[2,i,j] .* vals[i,j,2] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[3,2] == sum(vals[3,i,j] .* vals[i,j,2] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[1,3] == sum(vals[1,i,j] .* vals[i,j,3] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[2,3] == sum(vals[2,i,j] .* vals[i,j,3] for i in 1:3 for j in 1:3)
@test (t1 ⋅² t1)[3,3] == sum(vals[3,i,j] .* vals[i,j,3] for i in 1:3 for j in 1:3)

# a_il = b_ij*c_jl
v1 = [1 2 3
      2 4 0
      3 0 5];
v2 = [1 1 2
      1 2 5
      2 5 3];
t1 = TensorValue(v1)
t2 = TensorValue(v2)
@test (t1 ⋅ t2)[1,1] == sum(v1[1,j] .* v2[j,1] for j in 1:3)
@test (t1 ⋅ t2)[2,1] == sum(v1[2,j] .* v2[j,1] for j in 1:3)
@test (t1 ⋅ t2)[3,1] == sum(v1[3,j] .* v2[j,1] for j in 1:3)
@test (t1 ⋅ t2)[1,2] == sum(v1[1,j] .* v2[j,2] for j in 1:3)
@test (t1 ⋅ t2)[2,2] == sum(v1[2,j] .* v2[j,2] for j in 1:3)
@test (t1 ⋅ t2)[3,2] == sum(v1[3,j] .* v2[j,2] for j in 1:3)
@test (t1 ⋅ t2)[1,3] == sum(v1[1,j] .* v2[j,3] for j in 1:3)
@test (t1 ⋅ t2)[2,3] == sum(v1[2,j] .* v2[j,3] for j in 1:3)
@test (t1 ⋅ t2)[3,3] == sum(v1[3,j] .* v2[j,3] for j in 1:3)

# Complex

a = 1.0 + 3.0*im
b = 4.0 - 3.0*im
@test outer(a,b) == a*b
@test inner(a,b) == a*b

@test real(VectorValue(1+1im)) == VectorValue(1)
@test real(VectorValue(1+1im, 1+1im)) == VectorValue(1, 1)
@test real(VectorValue(1+1im, 1+1im, 1+1im)) == VectorValue(1, 1, 1)

@test imag(VectorValue(1+1im)) == VectorValue(1)
@test imag(VectorValue(1+1im, 1+1im)) == VectorValue(1, 1)
@test imag(VectorValue(1+1im, 1+1im, 1+1im)) == VectorValue(1, 1, 1)

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


end # module OperationsTests
