module MappingInterfacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Inference

using BenchmarkTools
using FillArrays

# Mapping Interfaces

a = [3,2]
b = [2,1]
test_mapping(+,(a,b),a+b)

bm = Broadcasting(+)
cache = return_cache(bm,a,b)
# @btime evaluate!($cache,$bm,$a,$b)

# m = rand(2,2)
# test_mapping(m,(a,b),m)
# testitem(m,a,b)

# m = rand(2,2)
# test_mapping(m,(a,b),m)
# testitem(m,a,b)

# cs = map(op -> return_cache(op,a,b),(+,m))
# cs = return_caches((+,m),a,b)
# evaluate!(cs,(+,m),a,b) == (a+b,m)
# evaluate((+,m),a,b) == (a+b,m)

# @test map(op -> return_type(op,a,b),(+,m)) == (Array{Int64,1}, Array{Float64,2})

z = evaluate(+,a,b)
c = return_cache(+,a,b)
evaluate!(c,+,a,b) == a+b
typeof(z) == typeof(a+b)
return_type(+,a,b)
testitem(+,a,b)


f = Broadcasting(+)
a = rand(3,2)
b = 3
c = a .+ b
Arrays.test_mapping(f,(a,b),c)

k = Broadcasting(-)
test_mapping(k,(1,),-1)
test_mapping(k,([1,2],),[-1,-2])
test_mapping(k,(1,2),-1)
test_mapping(k,(1.0,2),-1.0)
test_mapping(k,(1,2.0),-1.0)
test_mapping(k,([1,2],2),[-1,0])
test_mapping(k,(2,[1,2]),[1,0])
test_mapping(k,([3,4],[1,2]),[2,2])

f = Broadcasting(⋅)
a = fill(TensorValue(2,0,0,0,2,0,0,0,2),2)
b = VectorValue(1,2,3)
c = zeros(VectorValue{3,Int},2)
broadcast!(⋅,c,a,b)
test_mapping(f,(a,b),c)

cache = return_cache(f,a,b)
# @btime evaluate!($cache,$f,$a,$b)

# x = [1,2]
x = rand(3,3)

fa(x) = 2*x
test_mapping(fa,(x,),2*x)

fb(x) = sqrt.(x)
test_mapping(fb,(x,),sqrt.(x))

op = Broadcasting(*)
cache = return_cache(op,2,x)
# @btime evaluate!($cache,$op,$2,$x)
test_mapping(Broadcasting(*),(2,x),2*x)

fab = Operation(fa)(fb)
test_mapping(fab,(x,),2*(sqrt.(x)))

bm = Broadcasting(*)
cache = return_cache(bm,x,2)
# @btime evaluate!($cache,$bm,$x,$2)

bs = Broadcasting(sqrt)
cache = return_cache(bs,x)
# @btime evaluate!($cache,$bs,$x)

h = Operation(bs)(bm)
cache = return_cache(h,x,2)
# @btime evaluate!($cache,$h,$x,$2)


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


# k = Broadcasting(-)
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
