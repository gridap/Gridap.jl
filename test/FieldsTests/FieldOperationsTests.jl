module FieldOperationsTests

using Gridap.Arrays
using Gridap.Fields
using Gridap.TensorValues
using FillArrays

# using BenchmarkTools
using Test

# Operations

d = 2
p = Point(1.0,2.0)

np = 3
x = fill(p,np)

v = VectorValue{d}(1.0,1.0)
f = MockField(v)
fp = evaluate(f,p)

df = f+f
∇fp = 2*evaluate(gradient(f),p)
test_field(df,p,2*fp,grad=2.0*∇fp)

c = return_cache(df,p)
# @btime evaluate!($c,$df,$p)

∇df = ∇(df)
c = return_cache(∇df,p)
# @btime evaluate!($c,$∇df,$p)

df = f-f
test_field(df,p,fp-fp,grad=0.0*∇fp)

df = ConstantField(2.0)*f
test_field(df,p,fp*2.0,grad=2.0*∇fp)

c = return_cache(df,p)
# @btime evaluate!($c,$df,$p)

q(x) = 2*x
∇q = gradient(q)

f = GenericField(q)
fp = evaluate(f,p)
∇f = ∇(f)
∇fp = evaluate(∇f,p)
@test ∇fp == TensorValue(2.0, 0.0, 0.0, 2.0)
test_field(f,p,fp,grad=∇fp)

c = return_cache(f,p)
# @btime evaluate!($c,$f,$p)

df = f+f
test_field(df,p,2*fp,grad=2.0*∇fp)

c = return_cache(df,p)
# @btime evaluate!($c,$df,$p)

df = f-f
test_field(df,p,fp-fp,grad=0.0*∇fp)

df = ConstantField(2.0)*f
evaluate(df,p)
∇df = ∇(df)
∇dfp = 2*∇fp
evaluate(∇df,p)
test_field(df,p,fp*2.0,grad=∇dfp)

c = return_cache(df,p)
# @btime evaluate!($c,$df,$p)
c = return_cache(∇df,p)
# @btime evaluate!($c,$∇df,$p)


df = f⋅f
∇df = ∇(df)
dfp = (∇fp⋅fp)*2
test_field(df,p,fp⋅fp,grad=dfp)

c = return_cache(df,p)
# @btime evaluate!($c,$df,$p)
c = return_cache(∇df,p)
# @btime evaluate!($c,$∇df,$p)


bdf = df
evaluate(bdf,x)
∇bdf = gradient(bdf)
evaluate(∇bdf,x)
bdfp = fill(fp⋅fp,np)
∇bdfp = fill(dfp,np)
test_field(bdf,x,bdfp,grad=∇bdfp)

c = return_cache(bdf,x)
# @btime evaluate!($c,$bdf,$x)
c = return_cache(∇bdf,x)
# @btime evaluate!($c,$∇bdf,$x)

# Composition

q(x) = VectorValue(x[1]^2,x[2]^2)

f = GenericField(q)
bf = f
∇bf = ∇(bf)


evaluate(bf,x)
bdfx = evaluate(bf,evaluate(bf,x))
bdf = evaluate(Operation(bf),bf)
@test evaluate(bdf,x) == bdfx
test_field(bdf,x,bdfx)#,grad=∇bdfp)

bdfx = evaluate(bf+bf,x)
bdf = evaluate(Operation(+),bf,bf)
@test evaluate(bdf,x) == bdfx
test_field(bdf,x,bdfx)#,grad=∇bdfp)

bdfp = evaluate(bf⋅bf,p)
bdfx = evaluate(bf⋅bf,x)
bdf = bf⋅bf
@test evaluate(bdf,x) == bdfx
test_field(bdf,x,bdfx)#,grad=∇bdfp)

evaluate(bf,x)
evaluate(∇bf,x)
evaluate(bdf,x)
evaluate(∇bdf,x)
∇bdfx = evaluate(∇bf,x).⋅evaluate(bf,x)*2
evaluate(∇bdf,x) == ∇bdfx
test_field(bdf,x,bdfx,grad=∇bdfx)

# Composition

bdf = bf∘bf
∇bdf = ∇(bdf)
bdfx = evaluate(bdf,x)
∇bdfx = evaluate(∇bf,evaluate(bf,x)).⋅evaluate(∇bf,x)
test_field(bdf,x,bdfx,grad=∇bdfx)

c = return_cache(∇bdf,x)
# @btime evaluate!(c,∇bdf,x)

end # module
