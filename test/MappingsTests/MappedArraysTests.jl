module MappedArraysTests

using Test
using Gridap.Arrays
using Gridap.Mappings
using Gridap.Mappings: evaluate
using Gridap.Mappings: test_mapped_array
using FillArrays
using Gridap.Helpers

a = rand(3,2,4)
test_array(a,a)

a = rand(3,2,4)
a = CartesianIndices(a)
test_array(a,a)

a = rand(16,32)
a = CartesianIndices(a)
c = apply(-,a)
test_array(c,-a)

a = rand(12)
c = apply(-,a)
test_array(c,-a)

a = rand(12)
b = rand(12)
c = apply(-,a,b)
test_array(c,a.-b)

c = apply(Float64,-,a,b)
test_array(c,a.-b)

a = rand(0)
b = rand(0)
c = apply(-,a,b)
test_array(c,a.-b)

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
c = array_caches(a,b)
i = 1
v = getitems!(c,(a,b),i)
@test c == (nothing,nothing)
@test v == (a[i],b[i])

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
ai = testitem(a)
@test ai == a[1]
ai, bi = testitems(a,b)
@test ai == a[1]
@test bi == b[1]

a = fill(rand(Int,2,3),0)
b = fill(1,0)
ai = testitem(a)
@test ai == Array{Int,2}(undef,0,0)
ai, bi = testitems(a,b)
@test ai == Array{Int,2}(undef,0,0)
@test bi == zero(Int)

a = fill(+,10)
x = rand(10)
y = rand(10)
v = apply(a,x,y)
r = [(xi+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)
v = apply(Float64,a,x,y)
test_array(v,r)

p = 4
p1 = [1,2]
p2 = [2,1]
p3 = [4,3]
p4 = [6,1]

x = [p1,p2,p3,p4]

fa(x) = 2*x
fb(x) = sqrt.(x)

aa = Fill(fa,4)
r = apply(aa,x)
@test all([ r[i] ≈ 2*x[i] for i in 1:4])

bb = Fill(fb,4)
r = apply(bb,x)
@test all([ r[i] ≈ sqrt.(x[i]) for i in 1:4])

aaop = apply(operation,aa)
cm = apply(aaop,bb)
r = apply(cm,x)
@test all([ r[i] ≈ 2*(sqrt.(x[i])) for i in 1:4])

kk = cm[1]
ckk = return_cache(kk,p)
evaluate!(ckk,kk,p)


aop = apply(Operation(+),aa,bb)
apply(aa,x)+apply(bb,x)
apply(evaluate,aop,x)
@test apply(evaluate,aop,x) == apply(aa,x)+apply(bb,x)

# Broadcasting

a = fill(rand(2,3),12)
b = rand(12)
c = apply(Broadcasting(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),0)
b = rand(0)
c = apply(Broadcasting(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),12)
b = rand(12)
c = apply(Broadcasting(-),a,b)
d = apply(Broadcasting(+),a,c)
e = apply(Broadcasting(*),d,c)
test_array(e,[((ai.-bi).+ai).*(ai.-bi) for (ai,bi) in zip(a,b)])

a = Fill(Broadcasting(+),10)
x = [rand(2,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

a = Fill(Broadcasting(+),10)
x = [rand(mod(i-1,3)+1,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

# Operations

x = fill(4,3,3)

ax = Fill(x,4)
aa = Fill(fa,4)
bb = Fill(fb,4)

aop = apply(Operation(Broadcasting(+)),aa,bb)
aax = apply(evaluate,aa,ax)
bbx = apply(evaluate,bb,ax)
aopx = apply(evaluate,aop,ax)
@test aopx == aax+bbx

aop = apply(Operation(Broadcasting(*)),aa,bb)
aopx = apply(evaluate,aop,ax)
@test aopx[1] == aax[1].*bbx[1]

# Allocations

x = fill(4,3,3)

ax = Fill(x,4)
aa = Fill(Operation(fa),4)
bb = Fill(fb,4)
cm = apply(aa,bb)
r = apply(cm,ax)
@test all([ r[i] ≈ 2*(sqrt.(ax[i])) for i in 1:4])

nn = 2
an = Fill(nn,4)
ap = Fill(Broadcasting(*),4)
cm = apply(ap,ax,an)
@test all([cm[i] == nn*ax[i] for i in 1:4])

c_cm = Mappings.array_cache(cm)
@allocated Mappings.getindex!(c_cm,cm,1)
nalloc = 0
for i in length(cm)
  global nalloc
  nalloc += @allocated Mappings.getindex!(c_cm,cm,i)
end
@test nalloc == 0

as = Fill(Broadcasting(sqrt),4)
cs = apply(as,ax)
@test all([cs[i] == sqrt.(ax[i]) for i in 1:4])

c_cs = Mappings.array_cache(cs)
@allocated getindex!(c_cs,cs,1)
nalloc = 0
for i in length(cs)
  global nalloc
  nalloc += @allocated Mappings.getindex!(c_cs,cs,i)
end
@test nalloc == 0

asm = apply(operation,as)
ah = apply(asm,ap)
ch = apply(ah,ax,an)
@test all([ ch[i] ≈ sqrt.(nn*ax[i]) for i in 1:4])
c_ch = Mappings.array_cache(ch)
@allocated getindex!(c_ch,ch,1)
nalloc = 0
for i in length(ch)
  global nalloc
  nalloc += @allocated Mappings.getindex!(c_ch,ch,i)
end
@test nalloc == 0

end # module
