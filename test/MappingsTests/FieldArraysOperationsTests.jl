module FieldArraysOperationsTests

using Gridap.Arrays
using Gridap.Mappings
using Gridap.TensorValues
using FillArrays

using Test
using BenchmarkTools

# Check operations on results (OperationArray)

np = 4
d = 2
p = Point(1.0,2.0)
x = fill(p,np)
npp = np #(np,np)

v = VectorValue{d}(1.0,1.0)
f = MockField{d}(v)
fx = fill(v,np)

nf = 3
nff = nf #(nf,nf)

fa = fill(f,nf)

# Gradients

∇fa = BroadcastMapping(∇)(fa)
c = return_cache(∇fa,x)
# @santiagobadia : Not sure we can make it allocation-free
evaluate!(c,∇fa,x)
# @btime evaluate!($c,$∇fa,$x)

_∇fa = gradient.(fa)
c = return_cache(∇fa,x)
# @santiagobadia : Whereas this one is allocation-free...
# @btime evaluate!($c,$_∇fa,$x)

@test evaluate!(c,∇fa,x) == evaluate!(c,_∇fa,x)

# Transpose

c = return_cache(fa,p)
# @btime evaluate!($c,$fa,$p);

c = return_cache(fa,x)
# @btime evaluate!($c,$fa,$x);

tfa = transpose(fa)
c = return_cache(tfa,p)
# @btime evaluate!($c,$tfa,$p);
tr = evaluate!(c,tfa,p)
r = evaluate!(c,fa,p)
@test transpose(tr) == r

c = return_cache(tfa,x)
# @btime evaluate!($c,$tfa,$x);
tr = evaluate!(c,tfa,x)
r = evaluate!(c,fa,x)
@test size(tr) == (size(r,1),1,size(r,2))
@test reshape(permutedims(tr,[1,3,2]),size(r,1),size(r,2)) == r

# Broadcasted operations

op = +
b = BroadcastMapping(Operation(op))
bpfa = b(fa,fa)
@test bpfa == op(fa,fa)

fa isa AbstractArray{<:Field}
fa+fa
+(fa,fa)

c = return_cache(bpfa,x)
# @btime evaluate!($c,$bpfa,$x);

op = ⋅
b = BroadcastMapping(Operation(op))
bpfa = b(fa,fa)
@test bpfa == op.(fa,fa)

c = return_cache(bpfa,x)
# @btime evaluate!($c,$bpfa,$x);

# transpose

tfa = transpose(fa)
tr = evaluate(tfa,x)
r = evaluate(fa,x)

c = return_cache(tfa,x)
# @btime evaluate!($c,$tfa,$x);

# field + field array

op = ⋅
b = BroadcastMapping(Operation(op))
bpfa = b(f,fa)
@test evaluate(bpfa,x) == broadcast(⋅,evaluate(f,x),evaluate(fa,x))

c = return_cache(bpfa,x)
# @btime evaluate!(c,bpfa,x)

# column vector * row vector -> field array

op = ⋅
b = BroadcastMapping(Operation(op))
bpfa = b(fa,tfa)
evaluate(bpfa,x) == broadcast(⋅,evaluate(fa,x),evaluate(tfa,x))
@test all(bpfa .== broadcast(⋅,fa,tfa))

c = return_cache(bpfa,x)
# @btime evaluate!(c,bpfa,x)

# row vector * column vector -> field

dopa = tfa*fa
c = return_cache(dopa,x)
evaluate!(c,dopa,x)
@test evaluate(dopa,x)[1] == evaluate(tfa,p)⋅evaluate(fa,p)
# @btime evaluate!(c,dopa,x);

dopa = fa⋅fa
c = return_cache(dopa,x)
evaluate!(c,dopa,x)
@test evaluate(dopa,x)[1] == evaluate(tfa,p)⋅evaluate(fa,p)
# @btime evaluate!(c,dopa,x);

# linear combination

va = [1.0,2.0,3.0]
dopa = fa⋅va
@test evaluate(va⋅fa,x) == evaluate(fa⋅va,x)
c = return_cache(dopa,x)
# @btime evaluate!(c,dopa,x);

@test linear_combination(fa,va) == linear_combination(va,fa)
lc = linear_combination(fa,va)
fax = evaluate(fa,x)
lcx = evaluate(lc,x)
for i in 1:length(x)
  @test lcx[i] == fax[i,:]⋅va
end

c  =return_cache(lc,x)
# @btime evaluate!($c,$lc,$x);

vb = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
@test linear_combination(fa,vb) == linear_combination(vb,fa)
lc = linear_combination(fa,vb)
fax = evaluate(fa,x)
lcx = evaluate(lc,x)
for i in 1:length(x)
  @test lcx[i,:] == transpose(transpose(fax[i,:])*vb)
end

c  =return_cache(lc,x)
# @btime evaluate!($c,$lc,$x);

lc = linear_combination(fa,va)
fap = evaluate(fa,p)
lcx = evaluate(lc,p)
@test lcx == fap⋅va

c  =return_cache(lc,p)
# @btime evaluate!($c,$lc,$p);

lc = linear_combination(fa,vb)
fap = evaluate(fa,p)
lcp = evaluate(lc,p)
@test lcp == transpose(transpose(fap) * vb)

c  =return_cache(lc,p)
# @btime evaluate!($c,$lc,$p);

# composition

bm = BroadcastMapping(∘)
fc = bm(fa,f);

c = return_cache(fc,x)
# @btime evaluate!(c,fc,x)

end # module
