module FieldArraysTests

using Gridap.Fields
using Gridap.Arrays
using Gridap.TensorValues

using Test

# Testing the default interface for field arrays

function result(f,x) 
  T = return_type(testitem(f),testitem(x))
  r = zeros(T,size(x)...,size(f)...)
  for j in CartesianIndices(f)
    for i in CartesianIndices(x)
      r[i,j] = evaluate(f[j],x[i])
    end
  end
  r
end

p = Point(1.0,2.0)
np = 4
x = fill(p,np)
z = fill(p,0)

v = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
f = MockField.(v)

fp = v
∇fp = fill(zero(TensorValue{2,2,Float64}),length(v))
∇∇fp = fill(zero(ThirdOrderTensorValue{2,2,2,Float64,6}),length(v))
test_field_array(f,p,fp)
test_field_array(f,p,fp,grad=∇fp)
test_field_array(f,p,fp,grad=∇fp,gradgrad=∇∇fp)

test_field_array(f,x,result(f,x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x),gradgrad=result(∇∇.(f),x))

test_field_array(f,z,result(f,z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z),gradgrad=result(∇∇.(f),z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = Broadcasting(∇)(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)
#
#∇∇f = Broadcasting(∇∇)(f)
#c = return_cache(∇∇f,p)
#@btime evaluate!($c,$∇∇f,$p)
#c = return_cache(∇∇f,x)
#@btime evaluate!($c,$∇∇f,$x)

# testing 0-length array of fields

v = VectorValue{2,Float64}[]
f = MockField.(v)

fp = v
∇fp = fill(zero(TensorValue{2,2,Float64}),length(v))
∇∇fp = fill(zero(ThirdOrderTensorValue{2,2,2,Float64,6}),length(v))
test_field_array(f,p,fp)
test_field_array(f,p,fp,grad=∇fp)
test_field_array(f,p,fp,grad=∇fp,gradgrad=∇∇fp)

test_field_array(f,x,result(f,x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x),gradgrad=result(∇∇.(f),x))

test_field_array(f,z,result(f,z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z),gradgrad=result(∇∇.(f),z))

# Broadcasting operations

afun(x) = x[1]+2
a = GenericField(afun)

v = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
b = MockField.(v)

f = Broadcasting(Operation(*))(a,b)

fp = Broadcasting(*)(a(p),v)
∇fp = Broadcasting(⊗)(∇(a)(p),v)

test_field_array(f,p,fp)
test_field_array(f,p,fp,grad=∇fp)
test_field_array(f,x,result(f,x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x))
test_field_array(f,z,result(f,z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = Broadcasting(∇)(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)


# transpose of an array of fields

avals = [1.0,0.5,1.2]
a = MockField.(avals)

f = transpose(a)

fp = transpose(evaluate(a,p))
∇fp = transpose(evaluate(Broadcasting(∇)(a),p))

test_field_array(f,p,fp)
test_field_array(f,p,fp,grad=∇fp)
test_field_array(f,x,result(f,x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x))
test_field_array(f,z,result(f,z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z))

#using BenchmarkTools
#
#c = return_cache(f,p)
#@btime evaluate!($c,$f,$p)
#c = return_cache(f,x)
#@btime evaluate!($c,$f,$x)
#
#∇f = Broadcasting(∇)(f)
#c = return_cache(∇f,p)
#@btime evaluate!($c,$∇f,$p)
#c = return_cache(∇f,x)
#@btime evaluate!($c,$∇f,$x)

# transpose(i_to_f)*i_to_vals

avals = [1.0,0.5,1.2,3.9]
a = MockField.(avals)
b = VectorValue{2,Float64}[(1,1),(4,2),(3,5),(1,2)]

f = linear_combination(b,a)
@test isa(f,Fields.LinearCombinationField)

fp = transpose(b)*avals
∇fp = zero(TensorValue{2,2,Float64,4})

test_field(f,p,fp)
test_field(f,p,fp,grad=∇fp)
test_field(f,x,f.(x))
test_field(f,x,f.(x),grad=∇(f).(x))
test_field(f,z,f.(z))
test_field(f,z,f.(z),grad=∇(f).(z))


using BenchmarkTools

@btime linear_combination($b,$a)
c = return_cache(f,p)
@btime evaluate!($c,$f,$p)
c = return_cache(f,x)
@btime evaluate!($c,$f,$x)

∇f = Broadcasting(∇)(f)
c = return_cache(∇f,p)
@btime evaluate!($c,$∇f,$p)
c = return_cache(∇f,x)
@btime evaluate!($c,$∇f,$x)

# ij_to_vals*j_to_f

avals =rand(4)
a = MockField.(avals)
b = zeros(VectorValue{2,Float64},4,3)

f = linear_combination(b,a)
@test isa(f,Fields.LinearCombinationFieldVector)

fp = transpose(b)*avals
∇fp = zeros(TensorValue{2,2,Float64,4},3)

test_field_array(f,p,fp)
test_field_array(f,p,fp,grad=∇fp)
test_field_array(f,x,result(f,x))
test_field_array(f,x,result(f,x),grad=result(∇.(f),x))
test_field_array(f,z,result(f,z))
test_field_array(f,z,result(f,z),grad=result(∇.(f),z))

using BenchmarkTools
@btime linear_combination($b,$a)

c = return_cache(f,p)
@btime evaluate!($c,$f,$p)
c = return_cache(f,x)
@btime evaluate!($c,$f,$x)

∇f = Broadcasting(∇)(f)
c = return_cache(∇f,p)
@btime evaluate!($c,$∇f,$p)
c = return_cache(∇f,x)
@btime evaluate!($c,$∇f,$x)



## Test MockField
#
#np = 4
#p = Point(1.0,2.0)
#x = fill(p,np)
#npp = (np,np)
#
## v = 3.0
#d = 2
#v = VectorValue{d}(1.0,1.0)
#f = MockField{d}(v)
#fx = fill(v,np)
#
#nf = 3
#nff = (nf,nf)
#
#fa, p = test_field_array(f,p,nf,grad=true)
#
#c = return_cache(fa,p)
## @btime evaluate!($c,$fa,$p)
#c = return_cache(f,p)
## @btime evaluate!($c,$f,$p)
#
#bfa, x = test_broadcast_field_array(f,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x);
#
#∇f = gradient(f)
#
#∇fa, p = test_field_array(∇f,p,nf)
#bfa, x = test_broadcast_field_array(∇f,p,nf,np)
#
#@test gradient.(fa) == ∇fa
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#
#bf = f
#
#fa, x = test_field_array(bf,p,nf,grad=true)
#bfa, x = test_broadcast_field_array(bf,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(fa,p)
## @btime evaluate!($c,$fa,$p)
#
##
#
#q(x) = 2*x
#∇q = gradient(q)
#
#qf = GenericField(q)
#
#af, x = test_field_array(qf,p,nf,grad=true)
#bfa, x = test_broadcast_field_array(qf,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#
#bqf = qf
#
#af, x = test_field_array(bqf,p,nf,grad=true)
#bfa, x = test_broadcast_field_array(bqf,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#
#
## GenericField (constant)
#
#v = 1.0
## v = VectorValue(4.0,3.0)
#vf = GenericField(v)
#
#af, x = test_field_array(vf,p,nf,hessian=true)
#bfa, x = test_broadcast_field_array(vf,p,nf,np,hessian=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#c = return_cache(vf,p)
## @btime evaluate!($c,$vf,$p)
#
#v = VectorValue(4.0,3.0)
#vf = GenericField(v)
#
#af, x = test_field_array(vf,p,nf,grad=true)
#bfa, x = test_broadcast_field_array(vf,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#
## ZeroField
#
#zvf = ZeroField(vf)
#
#af, x = test_field_array(zvf,p,nf,grad=true)
#bfa, x = test_broadcast_field_array(zvf,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#
#bzvf = zvf
#
#af, x = test_field_array(bzvf,p,nf,grad=true)
#bfa, x = test_broadcast_field_array(bzvf,p,nf,np,grad=true)
#
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#
## Type unstable array
## @santiagobadia : Anything to do here?
## We would need to define promotion rules, convert, etc.
## It seems not feasible but I think this case is not of
## practical interest
## On the other hand, we can always optimize these methods
## for our particular problem at hand
#
#bfa = [zvf, vf, qf, zvf, vf, qf, zvf, vf, qf, zvf, vf, qf]
#
#c = return_cache(bfa,p)
## @btime evaluate!($c,$bfa,$p)
#
#
#q(x) = 2*x
#h(x) = 4*x
#qf = GenericField(q)
#hf = GenericField(h)
#
#bqf = qf
#bhf = hf
#
#af = [qf, hf, hf, qf]
#
#c = return_cache(af,p)
## @btime evaluate!($c,$af,$p)
#
#bfa = [bqf, bhf, bhf, bqf]
#c = return_cache(bfa,x)
## @btime evaluate!($c,$bfa,$x)

end # module
