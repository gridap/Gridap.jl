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
∇∇fp = zero(ThirdOrderTensorValue{2,2,2,Float64,8})
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

# integration

fun(x) = 3*x[1]
f = GenericField(fun)
ϕfun(x) = 2*x
ϕ = GenericField(ϕfun)

w = ones(size(x))

i = integrate(f,x,w)
@test i == sum(f.(x).*w)

i = integrate(f,x,w,∇(ϕ))
@test i == sum(f.(x).*w.*meas.(∇(ϕ).(x)))

#using BenchmarkTools
#c = return_cache(integrate,f,x,w)
#@btime evaluate!($c,$integrate,$f,$x,$w)
#
#J = ∇(ϕ)
#c = return_cache(integrate,f,x,w,J)
#@btime evaluate!($c,$integrate,$f,$x,$w,$J)


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
∇∇fp = zero(ThirdOrderTensorValue{2,2,2,Float64,8})
test_field(f,p,fp)
test_field(f,p,fp,grad=∇fp)
test_field(f,p,fp,grad=∇fp,gradgrad=∇∇fp)

test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,x,f.(x),grad=∇(f).(x),gradgrad=∇∇(f).(x))

test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))
test_field(f,z,f.(z),grad=∇(f).(z),gradgrad=∇∇(f).(z))

v = [1,2,3]
f = ConstantField.(v)
a = lazy_map(evaluate,f,fill([Point(1,2),Point(3,4)],length(v)))
test_array(a,[[1,1],[2,2],[3,3]])

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

f = Operation(/)(a,b)
cp = afun(p) / bfun(p)
test_field(f,p,cp)

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

Tfun(x) = diagonal_tensor(VectorValue(1*x[1],2*x[2]))
bfun(x) = VectorValue(x[2],x[1])
Fields.gradient(::typeof(Tfun)) = x-> ThirdOrderTensorValue{2,2,2,Float64}(1,0,0,0,0,0,0,2)
a = GenericField(Tfun)
b = GenericField(bfun)

f = Operation(⋅)(a,b)
cp = Tfun(p)⋅bfun(p)
∇cp = ∇(Tfun)(p)⋅bfun(p) + ∇(bfun)(p)⋅transpose(Tfun(p))
test_field(f,p,cp)
test_field(f,p,cp,grad=∇cp)
test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))

afun(x) = x.+2
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

afun(x) = x.+2
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

vfun(x) = 2*x[1]+x[2]
v = GenericField(vfun)
vt = VoidFieldMap(true)(v)
vf = VoidFieldMap(false)(v)
test_field(vt,p,zero(v(p)))
test_field(vf,p,v(p))
test_field(vt,p,zero(v(p)),grad=zero(∇(v)(p)))
test_field(vf,p,v(p),grad=∇(v)(p))
test_field(vt,p,zero(v(p)),grad=zero(∇(v)(p)),gradgrad=zero(∇∇(v)(p)))
test_field(vf,p,v(p),grad=∇(v)(p),gradgrad=∇∇(v)(p))
test_field(vt,x,zero.(v(x)))
test_field(vf,x,v(x))
test_field(vt,x,zero.(v(x)),grad=zero.(∇(v)(x)))
test_field(vf,x,v(x),grad=∇(v)(x))
test_field(vt,x,zero.(v(x)),grad=zero.(∇(v)(x)),gradgrad=zero.(∇∇(v)(x)))
test_field(vf,x,v(x),grad=∇(v)(x),gradgrad=∇∇(v)(x))

# testing hessian rule for sum and product of two fields

afun(x) = x[1]^3 + x[2]^4
bfun(x) = sin(x[1])*cos(x[2])
cfun(x) = exp(x⋅x)

a = GenericField(afun)
b = GenericField(bfun)
c = GenericField(cfun)

f = Operation(+)(Operation(*)(a,b), c)
∇f = ∇(a)*b + ∇(b)*a + ∇(c)
cp = afun(p) * bfun(p) + cfun(p)
∇cp = ∇(afun)(p)*bfun(p) + ∇(bfun)(p)*afun(p) + ∇(cfun)(p)
∇∇cp = ∇∇(afun)(p) * bfun(p) + afun(p) * ∇∇(bfun)(p) + ∇(afun)(p)⊗∇(bfun)(p) + ∇(bfun)(p)⊗∇(afun)(p) + ∇∇(cfun)(p)
test_field(f,p,cp)
test_field(f,p,cp, grad=∇cp, gradgrad=∇∇cp)

test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x),gradgrad=∇∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z),gradgrad=∇∇(f).(z))

# this one checks by taking ∇ of ∇f to see if matches with rule for ∇∇(f)
test_field(∇f, p, ∇cp, grad=∇∇cp)
test_field(∇f, x, ∇(f).(x), grad=∇∇(f).(x))
test_field(∇f, z, ∇(f).(z), grad=∇∇(f).(z))

# Product rules
#
# `gradient(::OperationField{typeof(op)})` delegates to `product_rule`, so each
# method is exercised through a real field operation and checked against central
# finite differences -- an oracle that assumes no index convention, which is
# exactly what these rules get wrong when they are wrong.

# ∂_k of a field by central differences, derivative index first as Gridap writes it
function fd(f,x::Point{D},h=1e-6) where D
  ntuple(Val(D)) do k
    e = VectorValue(ntuple(i -> i==k ? h : 0.0, Val(D)))
    (evaluate(f,x+e) - evaluate(f,x-e))/(2*h)
  end
end

# the value indices are contiguous within a slice of fixed k, since k leads and
# MultiValues are column major
function grad_matches(f,x,tol=1e-5)
  g = evaluate(gradient(f),x)
  d = fd(f,x)
  D = length(d)
  all(abs(g.data[(c-1)*D+k] - d[k].data[c]) < tol for k in 1:D for c in 1:length(d[k]))
end

function test_product_rule(op,f1,f2,x)
  h = Operation(op)(f1,f2)
  @test evaluate(h,x) ≈ op(evaluate(f1,x),evaluate(f2,x))
  @test grad_matches(h,x)
end

xp = Point(0.31,0.22)
sf = GenericField(x -> 1.0 + x[1]^2 + x[2])
vf = GenericField(x -> VectorValue(1.0 + x[1], 2.0 - x[2]^2))
wf = GenericField(x -> VectorValue(x[2], 1.0 + x[1]*x[2]))
mf = GenericField(x -> TensorValue(1.0 + x[1], 2.0*x[2], x[1]*x[2], 3.0 - x[2]))
nf = GenericField(x -> TensorValue(x[2], 1.0 - x[1], 2.0 + x[1]^2, x[1]*x[2]))

test_product_rule(⋅,vf,mf,xp)
test_product_rule(⋅,mf,nf,xp)
test_product_rule(outer,vf,sf,xp)
test_product_rule(outer,sf,vf,xp)
test_product_rule(outer,mf,sf,xp)
test_product_rule(outer,sf,mf,xp)
test_product_rule(outer,wf,vf,xp)
test_product_rule(outer,vf,mf,xp)
test_product_rule(outer,mf,vf,xp)
test_product_rule(outer,mf,nf,xp)

# What the ⋅ rules are for: the chain rule for ∇∇ of u∘F differentiates a matrix
# product. Checked against ∂_i∂_j (u∘F)_c = ∂_iF_a ∂_jF_b (∂_a∂_b u_c).
Faff = AffineField(TensorValue(2.0,0.5,-0.3,3.0),VectorValue(0.1,0.2))
uv = GenericField(x -> VectorValue(x[1]^2 + x[2]^3, 2.0*x[1]*x[2]))
Huf = evaluate(Broadcasting(∇∇)(Operation(uv)(Faff)),[xp])[1]
@test Huf isa ThirdOrderTensorValue{2,2,2}
Jaff = evaluate(gradient(Faff),xp)
Hu = evaluate(∇∇(uv),evaluate(Faff,xp))
for c in 1:2, i in 1:2, j in 1:2
  @test Huf[i,j,c] ≈ sum(Jaff[i,a]*Hu[a,b,c]*Jaff[j,b] for a in 1:2, b in 1:2) atol=1e-10
end

# What the outer rules are for: scattering a basis over the components of a
# bigger value, which was not differentiable at all.
function test_outer_scatter(b,x)
  sb = Broadcasting(Operation(outer))(ConstantField(VectorValue(1.0,0.0)),b)
  @test size(evaluate(Broadcasting(∇)(sb),[x])) == (1,length(b))
  @test all(grad_matches(sb[j],x) for j in 1:length(b))
end

test_outer_scatter([GenericField(x -> 1.0 + x[1]), GenericField(x -> x[1]*x[2])],xp)
test_outer_scatter([GenericField(x -> VectorValue(1.0 + x[1], x[2]))],xp)

end # module
