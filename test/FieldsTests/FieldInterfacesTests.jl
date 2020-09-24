module FieldInterfacesTests

using Gridap.Mappings
using Gridap.TensorValues

using BenchmarkTools
using Test

# using Gridap.NewFields: MockField, MockBasis, OtherMockBasis
# using Gridap.Mappings

# Test MockField

np = 4
p = Point(1.0,2.0)
x = fill(p,np)

# v = 3.0
d = 2
v = VectorValue{d}(1.0,1.0)
f = MockField{d}(v)
fx = fill(v,np)

evaluate(f,p)
bf = BroadcastField(f)
fp = evaluate(f,p)
bfx = fill(fp,np)
evaluate(bf,x)
test_field(bf,x,bfx)

c = return_cache(bf,x)
@btime evaluate!(c,bf,x)

@test BroadcastField(bf) == bf

# test_field(f,p,evaluate(f,p),evaluate(∇(f),p))

# @btime evaluate(f,x)
# @btime evaluate(bf,x)

∇f = gradient(f)
∇fp = evaluate(∇f,p)
@test ∇fp == zero(TensorValue{2,2,Float64,4})
test_field(f,p,fp,grad=∇fp)

c = return_cache(∇f,p)
@btime evaluate!(c,∇f,p)

∇bf = gradient(bf)
@test ∇bf == BroadcastField(gradient(f))
∇bfx = fill(∇fp,4)
test_field(∇bf,x,∇bfx)
test_field(bf,x,bfx,grad=∇bfx)

c = return_cache(∇bf,x)
@btime evaluate!(c,∇bf,x)

# GenericField (function)

@test GenericField(bf) == bf == Field(bf)
@test return_type(GenericField,Float64) == GenericField{Float64}
@test return_type(GenericField,BroadcastField) == BroadcastField

q(x) = 2*x
∇q = gradient(q)


qf = GenericField(q)
qfp = q(p)

@btime evaluate!(nothing,$qf,$p)

@test return_cache(qf,p) == return_cache(q,p)
@test return_type(qf,p) == return_type(q,p)
@test typeof(gradient(qf)) <: Gradient

∇qf = gradient(qf)
∇qfp = ∇q(p)
evaluate(∇qf,p)
test_field(qf,p,qfp,grad=∇qfp)

c = return_cache(∇qf,p)
@btime evaluate!($c,$∇qf,$p)


bqf = BroadcastField(qf)
∇bqf = ∇(bqf)
bqfx = q.(x)
∇bqfx = ∇q.(x)
test_field(bqf,x,bqfx,grad=∇bqfx)

c = return_cache(bqf,p)
@btime evaluate!($c,$bqf,$p)

c = return_cache(∇bqf,p)
@btime evaluate!($c,$∇bqf,$p)

# GenericField (constant)

v = 1.0
# v = VectorValue(4.0,3.0)
vf = GenericField(v)
test_field(vf,p,v)

c = return_cache(vf,p)
@btime evaluate!($c,$vf,$p)

∇vf = ∇(vf)
# ∇vp = TensorValue(0.0,0.0,0.0,0.0)
∇vp = VectorValue(0.0,0.0)
test_field(vf,p,v,grad=∇vp)

c = return_cache(∇vf,p)
@btime evaluate!($c,$∇vf,$p)

∇∇vf = gradient(gradient(vf))
evaluate(∇∇vf,p)
∇∇vfp = TensorValue(0.0,0.0,0.0,0.0)
test_field(vf,p,v,grad=∇vp,hessian=∇∇vfp)

c = return_cache(∇∇vf,p)
@btime evaluate!($c,$∇∇vf,$p)

# ZeroField

zvf = ZeroField(vf)
zvfp = zero(return_type(vf,p))
∇zvf = ∇(zvf)
∇zvfp = zero(return_type(∇vf,p))
test_field(zvf,p,zvfp,grad=∇zvfp)

bzvf = BroadcastField(zvf)
bzvfx = fill(zvfp,np)
∇bzvf = ∇(bzvf)
∇bzvfx = fill(∇zvfp,np)
test_field(bzvf,x,bzvfx,grad=∇bzvfx)

end # module
