module MockFieldsTests

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
evaluate(bf,x) == bfx
test_field(bf,x,bfx)

# test_field(f,p,evaluate(f,p),evaluate(∇(f),p))

# @btime evaluate(f,x)
# @btime evaluate(bf,x)

∇f = gradient(f)
∇fp = evaluate(∇f,p)
@test ∇fp == zero(TensorValue{2,2,Float64,4})
test_field(f,p,fp,grad=∇fp)


∇bf = gradient(bf)
@test ∇bf == BroadcastField(gradient(f))
∇bfx = fill(∇fp,4)
test_field(∇bf,x,∇bfx)
test_field(bf,x,bfx,grad=∇bfx)

# GenericField (function)

@test GenericField(bf) == bf == Field(bf)
@test return_type(GenericField,Float64) == GenericField{Float64}
@test return_type(GenericField,BroadcastField) == BroadcastField

q(x) = 2*x
∇q = gradient(q)

qf = GenericField(q)
qfp = q(p)

@test return_cache(qf,p) == return_cache(q,p)
@test return_type(qf,p) == return_type(q,p)
@test typeof(gradient(qf)) <: Gradient

∇qf = gradient(qf)
∇qfp = ∇q(p)
evaluate(∇qf,p)
test_field(qf,p,qfp,grad=∇qfp)

bqf = BroadcastField(qf)
∇bqf = ∇(bqf)
bqfx = q.(x)
∇bqfx = ∇q.(x)
test_field(bqf,x,bqfx,grad=∇bqfx)

# GenericField (constant)

v = 1.0
# v = VectorValue(4.0,3.0)
vf = GenericField(v)
test_field(vf,p,v)

∇vf = ∇(vf)
# ∇vp = TensorValue(0.0,0.0,0.0,0.0)
∇vp = VectorValue(0.0,0.0)
test_field(vf,p,v,grad=∇vp)


∇∇vf = gradient(gradient(vf))
evaluate(∇∇vf,p)
∇∇vfp = TensorValue(0.0,0.0,0.0,0.0)
test_field(vf,p,v,grad=∇vp,hessian=∇∇vfp)

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

# @btime evaluate(∇f,p)1
# xl = fill(p,1000)
# cache = return_cache(∇bf,xl)
# @btime evaluate!(cache,∇bf,xl)

# xl = fill(p,10000)
# @btime ∇q.(xl)
# c = return_cache(∇bqf,xl)
# @btime evaluate!(c,∇bqf,xl)





# Te = map(Mappings.numbertype,(p,))
# Mappings.return_type(f,Te...)
# Mappings.testitem(Te...)
# cache = return_cache(f,Te...)
# Mappings.testitem!(cache,f,Te...)
# evaluate!(cache,f,Te...)
# Te


# ∇fx = fill(VectorValue(v,0.0),np)
# ∇f = gradient(f)

# test_field(∇f,x,∇fx)
# test_field(f,x,fx,grad=∇fx)

# ndof = 8
# b = MockBasis(d,v,ndof)
# bx = fill(v,np,ndof)
# ∇bx = fill(VectorValue(v,0.0),np,ndof)
# test_field(b,x,bx,grad=∇bx)

# b = OtherMockBasis(d,ndof)
# bx = fill(2*p,np,ndof)
# ∇bx = fill(TensorValue(2.0,0.0,0.0,2.0),np,ndof)
# test_field(b,x,bx,grad=∇bx)

end # module
