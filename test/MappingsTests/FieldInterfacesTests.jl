module FieldInterfacesTests

using Gridap.Mappings
using Gridap.TensorValues
using Gridap.Arrays: CachedArray

using BenchmarkTools
using Test

np = 4
p = Point(1.0,2.0)
x = fill(p,np)

# v = 3.0
d = 2
v = VectorValue{d}(1.0,1.0)
f = MockField{d}(v)
fp = evaluate(f,p)
test_field(f,p,fp)

fx = fill(fp,np)
test_field(f,x,fx)

c = return_cache(f,x)
@btime evaluate!(c,f,x)

∇f = gradient(f)
∇fp = evaluate(∇f,p)
@test ∇fp == zero(TensorValue{2,2,Float64,4})
test_field(f,p,fp,grad=∇fp)

c = return_cache(∇f,p)
@btime evaluate!(c,∇f,p)

# GenericField (function)

@test return_type(GenericField,Float64) == GenericField{Float64}

q(x) = 2*x
∇q = gradient(q)


qf = GenericField(q)
qfp = q(p)

@btime evaluate!(nothing,$qf,$p)

@test return_cache(qf,p) == return_cache(q,p)
@test return_type(qf,p) == return_type(q,p)
@test typeof(gradient(qf)) <: FieldGradient

∇qf = gradient(qf)
∇qfp = ∇q(p)
evaluate(∇qf,p)
test_field(qf,p,qfp,grad=∇qfp)

c = return_cache(∇qf,p)
@btime evaluate!($c,$∇qf,$p)

qfx = q.(x)
∇qfx = ∇q.(x)
test_field(qf,x,qfx,grad=∇qfx)

c = return_cache(qf,p)
@btime evaluate!($c,$qf,$p)

c = return_cache(∇qf,p)
@btime evaluate!($c,$∇qf,$p)

c = return_cache(qf,x)
@btime evaluate!($c,$qf,$x)

c = return_cache(∇qf,x)
@btime evaluate!($c,$∇qf,$x)

# GenericField (constant)

v = 1.0
# v = VectorValue(4.0,3.0)
vf = GenericField(v)
test_field(vf,p,v)

vx = fill(v,np)
evaluate(vf,x)

test_field(vf,x,vx)

c = return_cache(vf,p)
@btime evaluate!($c,$vf,$p)

c = return_cache(vf,x)
@btime evaluate!($c,$vf,$x)


∇vf = ∇(vf)
∇vf.object
∇vp = VectorValue(0.0,0.0)
test_field(vf,p,v,grad=∇vp)
∇vfx = fill(∇vp,np)
test_field(vf,x,vx,grad=∇vfx)

c = return_cache(∇vf,x)
@btime evaluate!($c,$∇vf,$x)

c = return_cache(∇vf,x)
@btime evaluate!($c,$∇vf,$x)

∇∇vf = gradient(gradient(vf))
evaluate(∇∇vf,p)
∇∇vfp = TensorValue(0.0,0.0,0.0,0.0)
test_field(vf,p,v,grad=∇vp,hessian=∇∇vfp)

∇∇vfx = fill(∇∇vfp,np)
test_field(vf,x,vx,grad=∇vfx,hessian=∇∇vfx)

c = return_cache(∇∇vf,p)
@btime evaluate!($c,$∇∇vf,$p)

# ZeroField

zvf = ZeroField(vf)
zvfp = zero(return_type(vf,p))
∇zvf = ∇(zvf)
∇zvfp = zero(return_type(∇vf,p))
test_field(zvf,p,zvfp,grad=∇zvfp)

zvfx = fill(zvfp,np)
∇zvfx = fill(∇zvfp,np)

c = return_cache(∇zvf,x)

evaluate(zvf,x) == zvfx
test_field(zvf,x,zvfx,grad=∇zvfx)

end # module
