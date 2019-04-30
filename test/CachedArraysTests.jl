module CachedArraysTests

using Test
using Numa.CachedArrays

x = rand(10,5)

y = CachedArray(x)

@test y == x

z = y + x

a = y.array

setsize!(y,(5,3))

b = x[1:5,1:3]

@test y.array === a

@test size(y) == size(b)

@test y == b

setsize!(y,(15,30))

@test ! ( y.array === a )

@test size(y) == size(y.array)

y = CachedArray(Float64,2)

@test isa(y, AbstractArray{Float64,2})

@test isa(y, CachedArray{Float64,2,Array{Float64,2}})

@test size(y) == (0,0)

setsize!(y,(15,30))

@test size(y) == (15,30)

end # module CachedArraysTests
