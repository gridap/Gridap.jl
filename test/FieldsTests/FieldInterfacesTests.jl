module FieldInterfacesTests

using Gridap.Fields
using Gridap.TensorValues
using Gridap.Arrays

using Test

# Testing the default interface at a single point

p = Point(1.0,2.0)

# v = 3.0
v = VectorValue(1.0,1.0)
f = MockField(v)
fp = v
∇fp = zero(TensorValue{2,2,Float64})
∇∇fp = zero(ThirdOrderTensorValue{2,2,2,Float64,6})
test_field(f,p,fp)
test_field(f,p,fp,grad=∇fp)
test_field(f,p,fp,grad=∇fp,gradgrad=∇∇fp)

# Testing the default interface at a vector of points

np = 4
x = fill(p,np)
z = fill(p,0)

test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,x,f.(x),grad=∇(f).(x),gradgrad=∇∇(f).(x))

test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))
test_field(f,z,f.(z),grad=∇(f).(z),gradgrad=∇∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)
#
#∇∇f = ∇∇(f)
#c = return_cache(∇∇f,p)
#@btime evaluate!($c,$∇∇f,$p)
#c = return_cache(∇∇f,x)
#@btime evaluate!($c,$∇∇f,$x)

# Test field as collection

@test f === f[1]
_f, = f
@test f === _f
@test length(f) == 1
@test size(f) == ()
@test eltype(f) == typeof(f)

# GenericField (function)

q(x) = 2*x[1]

f = GenericField(q)

test_field(f,p,q(p))
test_field(f,p,q(p),grad=∇(q)(p))
test_field(f,p,q(p),grad=∇(q)(p),gradgrad=∇∇(q)(p))

test_field(f,x,q.(x))
test_field(f,x,q.(x),grad=∇(q).(x))
test_field(f,x,q.(x),grad=∇(q).(x),gradgrad=∇∇(q).(x))

test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))
test_field(f,z,f.(z),grad=∇(f).(z),gradgrad=∇∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)
#
#∇∇f = ∇∇(f)
#c = return_cache(∇∇f,p)
#@btime evaluate!($c,$∇∇f,$p)
#c = return_cache(∇∇f,x)
#@btime evaluate!($c,$∇∇f,$x)

# ZeroField

f = zero(f)
@test isa(f,ZeroField)

test_field(f,p,0*q(p))
test_field(f,p,0*q(p),grad=0*∇(q)(p))
test_field(f,p,0*q(p),grad=0*∇(q)(p),gradgrad=0*∇∇(q)(p))

test_field(f,x,0*q.(x))
test_field(f,x,0*q.(x),grad=0*∇(q).(x))
test_field(f,x,0*q.(x),grad=0*∇(q).(x),gradgrad=0*∇∇(q).(x))

test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))
test_field(f,z,f.(z),grad=∇(f).(z),gradgrad=∇∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)
#
#∇∇f = ∇∇(f)
#c = return_cache(∇∇f,p)
#@btime evaluate!($c,$∇∇f,$p)
#c = return_cache(∇∇f,x)
#@btime evaluate!($c,$∇∇f,$x)

# GenericField (function with more challenging domain)

h(x) = sqrt(x[1]-one(x[1]))
Arrays.testargs(::typeof(h),x) = (Point(map(one,x.data)),)

f = GenericField(h)
@test return_value(f,Point(0,0)) == 0.0
@test return_value(f,fill(Point(0,0),3)) == fill(0.0,3)

return_value(∇(f),Point(0,0))
return_value(∇(f),fill(Point(0,0),3))

return_value(∇∇(f),Point(0,0))
return_value(∇∇(f),fill(Point(0,0),3))

# ConstantField

v = VectorValue(1.0,1.0)
f = ConstantField(v)

fp = v
∇fp = zero(TensorValue{2,2,Float64})
∇∇fp = zero(ThirdOrderTensorValue{2,2,2,Float64,6})
test_field(f,p,fp)
test_field(f,p,fp,grad=∇fp)
test_field(f,p,fp,grad=∇fp,gradgrad=∇∇fp)

test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,x,f.(x),grad=∇(f).(x),gradgrad=∇∇(f).(x))

test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))
test_field(f,z,f.(z),grad=∇(f).(z),gradgrad=∇∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)
#
#∇∇f = ∇∇(f)
#c = return_cache(∇∇f,p)
#@btime evaluate!($c,$∇∇f,$p)
#c = return_cache(∇∇f,x)
#@btime evaluate!($c,$∇∇f,$x)

# Operations

afun(x) = x[1]+2
bfun(x) = sin(x[1])*cos(x[2])

a = GenericField(afun)
b = GenericField(bfun)

f = Operation(*)(a,b)
cp = afun(p) * bfun(p)
∇cp = ∇(afun)(p) * bfun(p) + afun(p) * ∇(bfun)(p)
test_field(f,p,cp)
test_field(f,p,cp,grad=∇cp)
test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)

afun(x) = x+2
bfun(x) = 2*x

a = GenericField(afun)
b = GenericField(bfun)

f = Operation(⋅)(a,b)
cp = afun(p)⋅bfun(p)
∇cp = ∇(afun)(p)⋅bfun(p) + ∇(bfun)(p)⋅afun(p)
test_field(f,p,cp)
test_field(f,p,cp,grad=∇cp)
test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)

afun(x) = x+2
bfun(x) = 2*x

a = GenericField(afun)
b = GenericField(bfun)

f = Operation(+)(a,b)
cp = afun(p)+bfun(p)
∇cp = ∇(afun)(p) + ∇(bfun)(p)
test_field(f,p,cp)
test_field(f,p,cp,grad=∇cp)
test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)

# Composition

mfun(g) = 3*g[1]
gfun(x) = 2*x
ffun(x) = mfun(gfun(x))

m = GenericField(mfun)
g = GenericField(gfun)

f = m∘g
fp = m(g(p))
∇fp = ∇(g)(p)⋅∇(m)(g(p))
test_field(f,p,fp)
test_field(f,p,fp,grad=∇fp)
test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = ∇(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)

#kk
#
#
#
#
#qfp = q(p)
#
## @btime evaluate!(nothing,$qf,$p)
#
#@test return_cache(qf,p) == return_cache(q,p)
#@test return_type(qf,p) == return_type(q,p)
#@test typeof(gradient(qf)) <: FieldGradient
#
#∇qf = gradient(qf)
#∇qfp = ∇q(p)
#evaluate(∇qf,p)
#test_field(qf,p,qfp,grad=∇qfp)
#
#c = return_cache(∇qf,p)
## @btime evaluate!($c,$∇qf,$p)
#
#qfx = q.(x)
#∇qfx = ∇q.(x)
#test_field(qf,x,qfx,grad=∇qfx)
#
#c = return_cache(qf,p)
## @btime evaluate!($c,$qf,$p)
#
#c = return_cache(∇qf,p)
## @btime evaluate!($c,$∇qf,$p)
#
#c = return_cache(qf,x)
## @btime evaluate!($c,$qf,$x)
#
#c = return_cache(∇qf,x)
## @btime evaluate!($c,$∇qf,$x)
#
## GenericField (constant)
#
#v = 1.0
## v = VectorValue(4.0,3.0)
#vf = GenericField(v)
#test_field(vf,p,v)
#
#vx = fill(v,np)
#evaluate(vf,x)
#
#test_field(vf,x,vx)
#
#c = return_cache(vf,p)
## @btime evaluate!($c,$vf,$p)
#
#c = return_cache(vf,x)
## @btime evaluate!($c,$vf,$x)
#
#
#∇vf = ∇(vf)
#∇vf.object
#∇vp = VectorValue(0.0,0.0)
#test_field(vf,p,v,grad=∇vp)
#∇vfx = fill(∇vp,np)
#test_field(vf,x,vx,grad=∇vfx)
#
#c = return_cache(∇vf,x)
## @btime evaluate!($c,$∇vf,$x)
#
#c = return_cache(∇vf,x)
## @btime evaluate!($c,$∇vf,$x)
#
#∇∇vf = gradient(gradient(vf))
#evaluate(∇∇vf,p)
#∇∇vfp = TensorValue(0.0,0.0,0.0,0.0)
#test_field(vf,p,v,grad=∇vp,hessian=∇∇vfp)
#
#∇∇vfx = fill(∇∇vfp,np)
#test_field(vf,x,vx,grad=∇vfx,hessian=∇∇vfx)
#
#c = return_cache(∇∇vf,p)
## @btime evaluate!($c,$∇∇vf,$p)
#
## ZeroField
#
#zvf = ZeroField(vf)
#zvfp = zero(return_type(vf,p))
#∇zvf = ∇(zvf)
#∇zvfp = zero(return_type(∇vf,p))
#test_field(zvf,p,zvfp,grad=∇zvfp)
#
#zvfx = fill(zvfp,np)
#∇zvfx = fill(∇zvfp,np)
#
#c = return_cache(∇zvf,x)
#
#evaluate(zvf,x) == zvfx
#test_field(zvf,x,zvfx,grad=∇zvfx)

end # module
