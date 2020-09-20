module MappingsTests

using Test
using Gridap.Arrays: CachedArray
using Gridap.Arrays
using Gridap.Mappings
using Gridap.TensorValues
using Gridap.Inference

using BenchmarkTools
using FillArrays

# Mapping Interfaces

a = [3,2]
b = [2,1]
test_mapping(+,(a,b),a+b)

m = rand(2,2)
test_mapping(m,(a,b),m)
testitem(m,a,b)

m = rand(2,2)
test_mapping(m,(a,b),m)
testitem(m,a,b)

cs = return_caches((+,m),a,b)
evaluate!(cs,(+,m),a,b) == (a+b,m)
evaluate((+,m),a,b) == (a+b,m)

return_types((+,m),a,b) == (Array{Int64,1}, Array{Float64,2})
Mappings.testitems(a,b) == (a,b)
Mappings._split(a,b,a,b) == (a,(b,a,b))
Mappings.return_types((+,m),a,b)
Mappings.return_types((+,/),1,1) == (Int,Float64)

z = evaluate(+,a,b)
c = return_cache(+,a,b)
evaluate!(c,+,a,b) == a+b
typeof(z) == typeof(a+b)
return_type(+,a,b)
testitem(+,a,b)


f = BroadcastMapping(+)
a = rand(3,2)
b = 3
c = a .+ b
Mappings.test_mapping(f,(a,b),c)

# Splatting issue
# cache = return_cache(f,a,b)
# evaluate!(cache,f,a,b)
# @test (@allocated evaluate!(cache,f,a,b)) == 0

k = BroadcastMapping(-)
test_mapping(k,(1,),-1)
test_mapping(k,([1,2],),[-1,-2])
test_mapping(k,(1,2),-1)
test_mapping(k,(1.0,2),-1.0)
test_mapping(k,(1,2.0),-1.0)
test_mapping(k,([1,2],2),[-1,0])
test_mapping(k,(2,[1,2]),[1,0])
test_mapping(k,([3,4],[1,2]),[2,2])

f = BroadcastMapping(⋅)
a = fill(TensorValue(2,0,0,0,2,0,0,0,2),2)
b = VectorValue(1,2,3)
c = zeros(VectorValue{3,Int},2)
broadcast!(⋅,c,a,b)
test_mapping(f,(a,b),c)

# x = [1,2]
x = rand(3,3)

fa(x) = 2*x
test_mapping(fa,(x,),2*x)

fb(x) = sqrt.(x)
test_mapping(fb,(x,),sqrt.(x))


fab = operation(fa,fb)
test_mapping(fab,(x,),2*(sqrt.(x)))

bm = BroadcastMapping(*)
c = return_cache(bm,x,2)
@allocated evaluate!(c,bm,x,2)
@test (@allocated evaluate!(c,bm,x,2)) == 0

bs = BroadcastMapping(sqrt)
d = return_cache(bs,x)
@allocated evaluate!(d,bs,x)
@test (@allocated evaluate!(d,bs,x)) == 0

h = operation(bs,bm)
e = return_cache(h,x,2)
@allocated evaluate!(e,h,x,2)
@test (@allocated evaluate!(e,h,x,2)) == 0

x = fill(4,3,3)

ax = Fill(x,4)
aa = Fill(fa,4)
bb = Fill(fb,4)
cm = apply(operation,aa,bb)
r = apply(cm,ax)
@test all([ r[i] ≈ 2*(sqrt.(ax[i])) for i in 1:4])

nn = 2
an = Fill(nn,4)
ap = Fill(BroadcastMapping(*),4)
cm = apply(ap,ax,an)
@test all([cm[i] == nn*ax[i] for i in 1:4])

c_cm = Mappings.array_cache(cm)
@allocated Mappings.getindex!(c_cm,cm,1)
nalloc = 0
for i in length(cm)
  global nalloc
  nalloc += @allocated Mappings.getindex!(c_cm,cm,i)
end
@test nalloc == 0

as = Fill(BroadcastMapping(sqrt),4)
cs = apply(as,ax)
@test all([cs[i] == sqrt.(ax[i]) for i in 1:4])

c_cs = Mappings.array_cache(cs)
@allocated getindex!(c_cs,cs,1)
nalloc = 0
for i in length(cs)
  global nalloc
  nalloc += @allocated Mappings.getindex!(c_cs,cs,i)
end
@test nalloc == 0

ah = apply(operation,as,ap)
ch = apply(ah,ax,an)
@test all([ ch[i] ≈ sqrt.(nn*ax[i]) for i in 1:4])
c_ch = Mappings.array_cache(ch)
@allocated getindex!(c_ch,ch,1)
nalloc = 0
for i in length(ch)
  global nalloc
  nalloc += @allocated Mappings.getindex!(c_ch,ch,i)
end
@test nalloc == 0

# k = -

# v = 3.0
# d = 2
# f = MockField{d}(v)
# g = composition(k,f)

# p = 4
# p1 = Point(1,2)
# p2 = Point(2,1)
# p3 = Point(4,3)
# p4 = Point(6,1)
# x = [p1,p2,p3,p4]

# fx = evaluate(f,x)
# ∇fx = evaluate(∇(f),x)
# gx = evaluate(k,fx)
# ∇gx = evaluate(k,∇fx)
# test_mapping(g,(x,),gx) #,grad=∇gx)

# # @santiagobadia : Create the gradient method for composition
# # using chain rule
# ∇g = composition(k,∇(f))
# test_mapping(∇g,(x,),∇gx)


# g = field_composition(k,f)
# @test g isa NewField
# gx = evaluate(g,x)
# cache = return_cache(g,x)
# test_mapping(g,(x,),gx) #grad=∇gx)

# fi = 3.0
# gi = 4.5
# d = 2
# f = MockField{d}(fi)
# g = MockField{d}(gi)

# # @santiagobadia : I like the tuple, k ∘ (f,g)
# h = composition(k,f,g)

# hx = evaluate(k,evaluate(f,x),evaluate(g,x))
# ∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
# test_mapping(h,(x,),hx) #,grad=∇hx)

# ∇h = composition(k,∇(f),∇(g))
# test_mapping(∇h,(x,),∇hx)

# fi = VectorValue(3.0,0.0)
# gi = VectorValue(4,5)
# d = 2
# f = MockField{d}(fi)
# g = ConstantField(gi)
# h = composition(k,f,g)

# hx = evaluate(k,evaluate(f,x),evaluate(g,x))
# ∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
# test_mapping(h,(x,),hx)

# ∇h = composition(k,∇(f),∇(g))
# test_mapping(∇h,(x,),∇hx)


# k = BroadcastMapping(-)
# fi = VectorValue(3.0,0.0)
# d = 2
# f = MockField{d}(fi)
# g = broadcast(ConstantField,x)
# h = composition(k,f,g)

# # @santiagobadia : A field can be paralelised wrt array of points
# # A mapping on the contrary can take multiple arguments...
# # We must create two evaluates, one for arrays of tuples and
# # one for one tuple
# # @santiagobadia : The previous version of gradient is wrong
# # @santiagobadia : We could also have derivative (transpose of gradient),
# # i.e., a Jacobian for f : Rn -> Rm

# hx = evaluate(k,evaluate(f,x),evaluate(g,x))
# test_mapping(h,(x,),hx) #,grad=∇hx)

# G = ∇(g)

# v = 3.0
# v = VectorValue(3.0,0.0)
# f = MockField{d}(v)

# f = ConstantField(v)

# fx = evaluate(f,x)
# test_field(f,(x,),fx)

end # module
