module FieldArraysOperationsTests

using Gridap.Arrays
using Gridap.Fields
using Gridap.TensorValues
using FillArrays

using Test
# using BenchmarkTools

np = 4
ncells = 10
d = 2
p = Point(1.0,2.0)
x = fill(p,np)
cell_to_x = Fill(x,ncells)


# ConstantField


cell_to_a = rand(ncells)
cell_to_f = lazy_map(ConstantField,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,map(i->fill(i,np),cell_to_a))

#print_op_tree(cell_to_fx)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)

# linear_combination -> Field

avals = [1.0,0.5,1.2,3.9]
a = MockField.(avals)
b = VectorValue{2,Float64}[(1,1),(4,2),(3,5),(1,2)]
f = linear_combination(b,a)

cell_to_a = Fill(a,ncells)
cell_to_b = fill(b,ncells)
cell_to_r = fill(f(x),ncells)

cell_to_f = lazy_map(linear_combination,cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
@test isa(cell_to_fx.maps.value,Fields.LinearCombinationMap)
test_array(cell_to_fx,cell_to_r)

cell_to_∇r = fill(∇(f)(x),ncells)
cell_to_∇f = lazy_map(∇,cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
@test isa(cell_to_∇fx.maps.value,Fields.LinearCombinationMap)
test_array(cell_to_∇fx,cell_to_∇r)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)

# linear_combination -> Field array

avals = [1.0,0.5,1.2,3.9]
a = MockField.(avals)
b = rand(length(avals),5)
f = linear_combination(b,a)

cell_to_a = Fill(a,ncells)
cell_to_b = fill(b,ncells)
cell_to_r = fill( evaluate(linear_combination(b,a),x),ncells)

cell_to_f = lazy_map(linear_combination,cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
@test isa(cell_to_fx.maps.value,Fields.LinearCombinationMap)
test_array(cell_to_fx,cell_to_r)

cell_to_∇r = fill(evaluate(Broadcasting(∇)(f),x),ncells)
cell_to_∇f = lazy_map(Broadcasting(∇),cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
@test isa(cell_to_∇fx.maps.value,Fields.LinearCombinationMap)
test_array(cell_to_∇fx,cell_to_∇r)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)


# transpose of an array of fields

avals = [1.0,0.5,1.2,3.9]
a = MockField.(avals)
b = rand(length(avals),5)
c = linear_combination(b,a)

cell_to_a = Fill(a,ncells)
cell_to_b = fill(b,ncells)
cell_to_c = lazy_map(linear_combination,cell_to_b,cell_to_a)

cell_to_f = lazy_map(transpose,cell_to_c)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
@test isa(cell_to_fx.maps.value,Fields.TransposeMap)
cell_to_r = fill(evaluate(transpose(c),x),ncells)
test_array(cell_to_fx,cell_to_r)

cell_to_∇f = lazy_map(Broadcasting(∇),cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
@test isa(cell_to_fx.maps.value,Fields.TransposeMap)
cell_to_∇r = fill(evaluate(transpose(∇.(c)),x),ncells)
test_array(cell_to_∇fx,cell_to_∇r)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)

# Composition between fields

mfun(g) = 3*g[1]
gfun(x) = 2*x
ffun(x) = mfun(gfun(x))

m = GenericField(mfun)
g = GenericField(gfun)
f = m∘g

cell_to_m = Fill(m,ncells)
cell_to_g = fill(g,ncells)
cell_to_f = lazy_map(∘,cell_to_m,cell_to_g)
cell_to_r = fill(f(x),ncells)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)

cell_to_∇f = lazy_map(∇,cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
cell_to_∇r = fill(evaluate(∇(f),x),ncells)
test_array(cell_to_∇fx,cell_to_∇r)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)

# Composition between field array and field

afun(x) = 2*x
bfun(x) = x[1] + 3
a = GenericField(afun)
b = fill(GenericField(bfun),4)
f = Broadcasting(∘)(b,a)

cell_to_a = fill(a,ncells)
cell_to_b = Fill(b,ncells)

cell_to_f = lazy_map(Broadcasting(∘),cell_to_b,cell_to_a)
cell_to_r = fill(evaluate(f,x),ncells)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)

cell_to_∇f = lazy_map(Broadcasting(∇),cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
cell_to_∇r = fill(evaluate(∇.(f),x),ncells)
test_array(cell_to_∇fx,cell_to_∇r)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)

# Operations

afun(x) = x[2]
bfun(x) = x[1] + 3
a = GenericField(afun)
b = GenericField(bfun)
f = Operation(+)(a,b)

cell_to_a = fill(a,ncells)
cell_to_b = Fill(b,ncells)
cell_to_r = fill(evaluate(f,x),ncells)
cell_to_∇r = fill(evaluate(∇(f),x),ncells)

cell_to_f = lazy_map(Operation(+),cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)
@test cell_to_fx.maps.value == Broadcasting(+)

cell_to_∇f = lazy_map(∇,cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
test_array(cell_to_∇fx,cell_to_∇r)
@test cell_to_fx.maps.value == Broadcasting(+)

T = GenericField{Nothing}
cell_to_f = lazy_map(Operation(+),T,cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)
@test cell_to_fx.maps.value == Broadcasting(+)

cell_to_∇f = lazy_map(∇,cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
test_array(cell_to_∇fx,cell_to_∇r)
@test cell_to_fx.maps.value == Broadcasting(+)

#using BenchmarkTools
#c = array_cache(cell_to_fx)
#@btime getindex!($c,$cell_to_fx,1)

# Operations with Broadcasting

afun(x) = x[2]
bfun(x) = x[1] + 3
a = GenericField(afun)
b = fill(GenericField(bfun),4)
f = Broadcasting(Operation(+))(a,b)

cell_to_a = fill(a,ncells)
cell_to_b = Fill(b,ncells)
cell_to_r = fill(evaluate(f,x),ncells)
cell_to_∇r = fill(evaluate(Broadcasting(∇)(f),x),ncells)

cell_to_f = lazy_map(Broadcasting(Operation(+)),cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)
#@test cell_to_fx.maps.value == Broadcasting(+)
@test cell_to_fx.maps.value == Fields.BroadcastingFieldOpMap(+)

cell_to_∇f = lazy_map(Broadcasting(∇),cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
test_array(cell_to_∇fx,cell_to_∇r)
#@test cell_to_∇fx.maps.value == Broadcasting(+)
@test cell_to_∇fx.maps.value == Fields.BroadcastingFieldOpMap(+)

T = GenericField{Nothing}
cell_to_f = lazy_map(Broadcasting(Operation(+)),T,cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)
#@test cell_to_fx.maps.value == Broadcasting(+)
@test cell_to_fx.maps.value == Fields.BroadcastingFieldOpMap(+)

# Operations with Broadcasting (product)

afun(x) = x[2]
bfun(x) = x[1] + 3
a = GenericField(afun)
b = fill(GenericField(bfun),4)
f = Broadcasting(Operation(*))(a,b)

cell_to_a = fill(a,ncells)
cell_to_b = Fill(b,ncells)
cell_to_r = fill(evaluate(f,x),ncells)
cell_to_∇r = fill(evaluate(Broadcasting(∇)(f),x),ncells)

cell_to_f = lazy_map(Broadcasting(Operation(*)),cell_to_b,cell_to_a)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
test_array(cell_to_fx,cell_to_r)

cell_to_∇f = lazy_map(Broadcasting(∇),cell_to_f)
cell_to_∇fx = lazy_map(evaluate,cell_to_∇f,cell_to_x)
test_array(cell_to_∇fx,cell_to_∇r)

# Integration

avals = [1.0,0.5,1.2,3.9]
a = MockField.(avals)
b = rand(length(avals),5)
d = linear_combination(b,a)
c = Broadcasting(Operation(*))(d,transpose(d))
w = rand(np)
fx = integrate(c,x,w)

cell_to_a = Fill(a,ncells)
cell_to_b = fill(b,ncells)
cell_to_w = Fill(w,ncells)
cell_to_d = lazy_map(linear_combination,cell_to_b,cell_to_a)
cell_to_dt = lazy_map(transpose,cell_to_d)
cell_to_c = lazy_map(Broadcasting(Operation(*)),cell_to_d,cell_to_dt)
cell_to_fx = lazy_map(integrate,cell_to_c,cell_to_x,cell_to_w)
cell_to_r = fill(fx,ncells)
test_array(cell_to_fx,cell_to_r)
@test isa(cell_to_fx.maps.value,IntegrationMap)

ϕ = GenericField(x->2*x)
fx = integrate(c,x,w,∇(ϕ))

cell_to_ϕ = fill(ϕ,ncells)
cell_to_J = lazy_map(∇,cell_to_ϕ)

cell_to_fx = lazy_map(integrate,cell_to_c,cell_to_x,cell_to_w,cell_to_J)
cell_to_r = fill(fx,ncells)
test_array(cell_to_fx,cell_to_r)
@test isa(cell_to_fx.maps.value,IntegrationMap)

# Reindex

avals = [1.0,0.5,1.2,3.9]
a = MockField.(avals)
b = VectorValue{2,Float64}[(1,1),(4,2),(3,5),(1,2)]
f = linear_combination(b,a)

cell_to_a = Fill(a,ncells)
cell_to_b = fill(b,ncells)
cell_to_c = lazy_map(linear_combination,cell_to_b,cell_to_a)

j_to_cell = [2,3,2,2,1]
j_to_c = lazy_map(Reindex(cell_to_c),j_to_cell)
j_to_x = Fill(x,length(j_to_cell))

j_to_cx = lazy_map(evaluate,j_to_c,j_to_x)
j_to_r = fill( f(x), length(j_to_cell) )
@test isa(j_to_cx.maps.value,Fields.LinearCombinationMap)
test_array(j_to_cx,j_to_r)

j_to_∇c = lazy_map(Broadcasting(∇),j_to_c)
j_to_∇cx = lazy_map(evaluate,j_to_∇c,j_to_x)
@test isa(j_to_∇cx.maps.value,Fields.LinearCombinationMap)
j_to_∇r = fill( ∇(f)(x), length(j_to_cell) )
test_array(j_to_∇cx,j_to_∇r)

# PosNegReindex

aval = 1.0
a = MockField(aval)
ipos_to_cell = [2,5,1]
cell_to_i = PosNegPartition(ipos_to_cell,ncells)
npos = length(cell_to_i.ipos_to_i)
nneg = length(cell_to_i.ineg_to_i)

ipos_to_a = fill(a,npos)
ineg_to_a = fill(0*a,nneg)

T = Union{eltype(ipos_to_a),eltype(ineg_to_a)}
cell_to_f = lazy_map(PosNegReindex(ipos_to_a,ineg_to_a),T,cell_to_i)
cell_to_fx = lazy_map(evaluate,cell_to_f,cell_to_x)
cell_to_r = [ zeros(np) for cell in 1:ncells]
cell_to_r[ipos_to_cell] = [ fill(aval,np) for cell in 1:npos]
#print_op_tree(cell_to_fx)

@test isa(cell_to_fx.maps.value,PosNegReindex)
test_array(cell_to_fx,cell_to_r)

#npp = np #(np,np)
#
#v = VectorValue{d}(1.0,1.0)
#f = MockField{d}(v)
#fx = fill(v,np)
#
#nf = 3
#nff = nf #(nf,nf)
#
#fa = fill(f,nf)
#
## Gradients
#
#∇fa = Broadcasting(∇)(fa)
#c = return_cache(∇fa,x)
## @santiagobadia : Not sure we can make it allocation-free
#evaluate!(c,∇fa,x)
## @btime evaluate!($c,$∇fa,$x)
#
#_∇fa = gradient.(fa)
#c = return_cache(∇fa,x)
## @santiagobadia : Whereas this one is allocation-free...
## @btime evaluate!($c,$_∇fa,$x)
#
#@test evaluate!(c,∇fa,x) == evaluate!(c,_∇fa,x)
#
## Transpose
#
#c = return_cache(fa,p)
## @btime evaluate!($c,$fa,$p);
#
#c = return_cache(fa,x)
## @btime evaluate!($c,$fa,$x);
#
#tfa = transpose(fa)
#c = return_cache(tfa,p)
## @btime evaluate!($c,$tfa,$p);
#tr = evaluate!(c,tfa,p)
#r = evaluate!(c,fa,p)
#@test transpose(tr) == r
#
#c = return_cache(tfa,x)
## @btime evaluate!($c,$tfa,$x);
#tr = evaluate!(c,tfa,x)
#r = evaluate!(c,fa,x)
#@test size(tr) == (size(r,1),1,size(r,2))
#@test reshape(permutedims(tr,[1,3,2]),size(r,1),size(r,2)) == r
#
## Broadcasted operations
#
#op = +
#b = Broadcasting(Operation(op))
#bpfa = b(fa,fa)
#@test bpfa == op(fa,fa)
#
#fa isa AbstractArray{<:Field}
#fa+fa
#+(fa,fa)
#
#c = return_cache(bpfa,x)
## @btime evaluate!($c,$bpfa,$x);
#
#op = ⋅
#b = Broadcasting(Operation(op))
#bpfa = b(fa,fa)
#@test bpfa == op.(fa,fa)
#
#c = return_cache(bpfa,x)
## @btime evaluate!($c,$bpfa,$x);
#
## transpose
#
#tfa = transpose(fa)
#tr = evaluate(tfa,x)
#r = evaluate(fa,x)
#
#c = return_cache(tfa,x)
## @btime evaluate!($c,$tfa,$x);
#
## field + field array
#
#op = ⋅
#b = Broadcasting(Operation(op))
#bpfa = b(f,fa)
#@test evaluate(bpfa,x) == broadcast(⋅,evaluate(f,x),evaluate(fa,x))
#
#c = return_cache(bpfa,x)
## @btime evaluate!(c,bpfa,x)
#
## column vector * row vector -> field array
#
#op = ⋅
#b = Broadcasting(Operation(op))
#bpfa = b(fa,tfa)
#evaluate(bpfa,x) == broadcast(⋅,evaluate(fa,x),evaluate(tfa,x))
#@test all(bpfa .== broadcast(⋅,fa,tfa))
#
#c = return_cache(bpfa,x)
## @btime evaluate!(c,bpfa,x)
#
## row vector * column vector -> field
#
#dopa = tfa*fa
#c = return_cache(dopa,x)
#evaluate!(c,dopa,x)
#@test evaluate(dopa,x)[1] == evaluate(tfa,p)⋅evaluate(fa,p)
## @btime evaluate!(c,dopa,x);
#
#dopa = fa⋅fa
#c = return_cache(dopa,x)
#evaluate!(c,dopa,x)
#@test evaluate(dopa,x)[1] == evaluate(tfa,p)⋅evaluate(fa,p)
## @btime evaluate!(c,dopa,x);
#
## linear combination
#
#va = [1.0,2.0,3.0]
#dopa = fa⋅va
#@test evaluate(va⋅fa,x) == evaluate(fa⋅va,x)
#c = return_cache(dopa,x)
## @btime evaluate!(c,dopa,x);
#
#@test linear_combination(fa,va) == linear_combination(va,fa)
#lc = linear_combination(fa,va)
#fax = evaluate(fa,x)
#lcx = evaluate(lc,x)
#for i in 1:length(x)
#  @test lcx[i] == fax[i,:]⋅va
#end
#
#c  =return_cache(lc,x)
## @btime evaluate!($c,$lc,$x);
#
#vb = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
#@test linear_combination(fa,vb) == linear_combination(vb,fa)
#lc = linear_combination(fa,vb)
#fax = evaluate(fa,x)
#lcx = evaluate(lc,x)
#for i in 1:length(x)
#  @test lcx[i,:] == transpose(transpose(fax[i,:])*vb)
#end
#
#c  =return_cache(lc,x)
## @btime evaluate!($c,$lc,$x);
#
#lc = linear_combination(fa,va)
#fap = evaluate(fa,p)
#lcx = evaluate(lc,p)
#@test lcx == fap⋅va
#
#c  =return_cache(lc,p)
## @btime evaluate!($c,$lc,$p);
#
#lc = linear_combination(fa,vb)
#fap = evaluate(fa,p)
#lcp = evaluate(lc,p)
#@test lcp == transpose(transpose(fap) * vb)
#
#c  =return_cache(lc,p)
## @btime evaluate!($c,$lc,$p);
#
## composition
#
#bm = Broadcasting(∘)
#fc = bm(fa,f);
#
#c = return_cache(fc,x)
## @btime evaluate!(c,fc,x)

end # module
