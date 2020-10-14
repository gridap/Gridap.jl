module LazyArraysTests

using Test
using Gridap.Arrays
using FillArrays
using Gridap.Helpers

a = rand(3,2,4)
test_array(a,a)

a = rand(3,2,4)
a = CartesianIndices(a)
test_array(a,a)

a = rand(16,32)
a = CartesianIndices(a)
c = lazy_map(-,a)
test_array(c,-a)

a = rand(12)
c = lazy_map(-,a)
test_array(c,-a)

a = rand(12)
b = rand(12)
c = lazy_map(-,a,b)
test_array(c,a.-b)

c = lazy_map(-,Float64,a,b)
test_array(c,a.-b)

a = rand(0)
b = rand(0)
c = lazy_map(-,a,b)
test_array(c,a.-b)

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
c = map(array_cache,(a,b))
i = 1
v = map((ci,ai) -> getindex!(ci,ai,i),c,(a,b))
# v = getitems!(c,(a,b),i)
@test c == (nothing,nothing)
@test v == (a[i],b[i])

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
ai = testitem(a)
@test ai == a[1]
ai, bi = map(testitem,(a,b))
@test ai == a[1]
@test bi == b[1]

a = fill(rand(Int,2,3),0)
b = fill(1,0)
ai = testitem(a)
@test ai == Array{Int,2}(undef,0,0)
ai, bi = map(testitem,(a,b))
@test ai == Array{Int,2}(undef,0,0)
@test bi == zero(Int)

a = fill(+,10)
x = rand(10)
y = rand(10)
v = lazy_map(+,x,y)
r = [(xi+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)
v = lazy_map(+,Float64,x,y)
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
r = lazy_map(fa,x)
@test all([ r[i] ≈ 2*x[i] for i in 1:4])

bb = Fill(fb,4)
r = lazy_map(fb,x)
@test all([ r[i] ≈ sqrt.(x[i]) for i in 1:4])

aaop = lazy_map(operation,aa)
cm = lazy_map(evaluate,aaop,bb)
r = lazy_map(evaluate,cm,x)
@test all([ r[i] ≈ 2*(sqrt.(x[i])) for i in 1:4])

kk = cm[1]
ckk = return_cache(kk,p)
evaluate!(ckk,kk,p)


aop = lazy_map(Operation(+),aa,bb)
lazy_map(evaluate,aa,x)+lazy_map(evaluate,bb,x)
lazy_map(evaluate,aop,x)
@test lazy_map(evaluate,aop,x) == lazy_map(evaluate,aa,x)+lazy_map(evaluate,bb,x)

# Broadcasting

a = fill(rand(2,3),12)
b = rand(12)
c = lazy_map(Broadcasting(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),0)
b = rand(0)
c = lazy_map(Broadcasting(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),12)
b = rand(12)
c = lazy_map(Broadcasting(-),a,b)
d = lazy_map(Broadcasting(+),a,c)
e = lazy_map(Broadcasting(*),d,c)
test_array(e,[((ai.-bi).+ai).*(ai.-bi) for (ai,bi) in zip(a,b)])

a = Fill(Broadcasting(+),10)
x = [rand(2,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = lazy_map(evaluate,a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

a = Fill(Broadcasting(+),10)
x = [rand(mod(i-1,3)+1,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = lazy_map(evaluate,a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

# Operations

x = fill(4,3,3)

ax = Fill(x,4)
aa = Fill(fa,4)
bb = Fill(fb,4)

aop = lazy_map(Operation(Broadcasting(+)),aa,bb)
aax = lazy_map(evaluate,aa,ax)
bbx = lazy_map(evaluate,bb,ax)
aopx = lazy_map(evaluate,aop,ax)
@test aopx == aax+bbx

aop = lazy_map(Operation(Broadcasting(*)),aa,bb)
aopx = lazy_map(evaluate,aop,ax)
@test aopx[1] == aax[1].*bbx[1]

# Allocations

x = fill(4,3,3)

ax = Fill(x,4)
aa = Fill(Operation(fa),4)
bb = Fill(fb,4)
cm = lazy_map(evaluate,aa,bb)
r = lazy_map(evaluate,cm,ax)
@test all([ r[i] ≈ 2*(sqrt.(ax[i])) for i in 1:4])

nn = 2
an = Fill(nn,4)
ap = Fill(Broadcasting(*),4)
cm = lazy_map(evaluate,ap,ax,an)
@test all([cm[i] == nn*ax[i] for i in 1:4])

c_cm = array_cache(cm)
@allocated getindex!(c_cm,cm,1)
nalloc = 0
for i in length(cm)
  global nalloc
  nalloc += @allocated getindex!(c_cm,cm,i)
end
# @test nalloc == 0

as = Fill(Broadcasting(sqrt),4)
cs = lazy_map(evaluate,as,ax)
@test all([cs[i] == sqrt.(ax[i]) for i in 1:4])

c_cs = array_cache(cs)
@allocated getindex!(c_cs,cs,1)
nalloc = 0
for i in length(cs)
  global nalloc
  nalloc += @allocated getindex!(c_cs,cs,i)
end
# @test nalloc == 0

asm = lazy_map(operation,as)
ah = lazy_map(evaluate,asm,ap)
ch = lazy_map(evaluate,ah,ax,an)
@test all([ ch[i] ≈ sqrt.(nn*ax[i]) for i in 1:4])
c_ch = array_cache(ch)
@allocated getindex!(c_ch,ch,1)
nalloc = 0
for i in length(ch)
  global nalloc
  nalloc += @allocated getindex!(c_ch,ch,i)
end
# @test nalloc == 0


# Empty arrays

g(x) = 1.0*x
a = lazy_map(g,Int[])
@test test_array(a,Float64[])

f(x) = sqrt(x-1)
@test_broken test_array(lazy_map(f,Int[]),Float64[])

# Fill optimizations

l = 10
v = 3
a = lazy_map(g,Fill(v,l))
test_array(a,Fill(g(v),l))
@test isa(a,Fill)

a = lazy_map(g,Float64,Fill(v,l))
test_array(a,Fill(g(v),l))
@test isa(a,Fill)

# array_cache

a = [rand(3) for i in 1:l]
b = [rand(3) for i in 1:l]
c = lazy_map(Broadcasting(+),a,b)
d = lazy_map(Broadcasting(*),a,c)

cache = array_cache(d)
@test getindex!(cache,d,1) === getindex!(cache,d,2)

i = 1
cache = array_cache(d,i)
@test getindex!(cache,d,i) === getindex!(cache,d,i)

ci = CartesianIndex((1))
cache = array_cache(d,ci)
@test getindex!(cache,d,ci) === getindex!(cache,d,ci)

# Shapes and indices

s = (3,4)
a = rand(s...)
b = rand(s...)
c = lazy_map(+,a,b)
@test size(c) == s
d = a.+b
test_array(c,d)
@test c[1] == d[1]
@test c[2,3] == d[2,3]
@test c[CartesianIndex(2,3)] == d[CartesianIndex(2,3)]
@test c[2,3,1,1,1] == d[2,3,1,1,1]
@test c[CartesianIndex(2,3),1,1] == d[CartesianIndex(2,3),1,1]

s = (3,4)
a = rand(s...)
b = rand(prod(s))
c = lazy_map(+,a,b)
d = map(+,a,b)
@test size(c) == size(d)
@test  c == d


end # module
