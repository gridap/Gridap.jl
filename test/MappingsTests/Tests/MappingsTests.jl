module MappingsTests

using Test
using Gridap.Arrays: CachedArray
using Gridap.Mappings
using Gridap.TensorValues

# Mapping Interfaces

a = [3,2]
b = [2,1]
test_mapping(+,(a,b),a+b)
testitem(+,a,b)

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

x = [1,2]

fa(x) = 2*x
test_mapping(fa,(x,),[2,4])

fb(x) = sqrt.(x)
test_mapping(fb,(x,),[1,sqrt(2)])

end # module
