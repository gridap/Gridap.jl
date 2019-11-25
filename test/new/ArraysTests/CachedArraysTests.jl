module CachedArraysTests

using Test
using Gridap.Arrays

x = rand(10,5)

y = CachedArray(x)

@test y == x

z = y + x

a = y.array

setsize!(y,(5,3))
y .= 0
@test y == y.array

setsize!(y,(Int32(5),Int32(3)))
y .= 0
@test y == y.array


setsize!(y,(10,5))
@test y == y.array
@test x === y.array

setsize!(y,(5,3))
@test y == y.array

b = x[1:5,1:3]

@test size(y) == size(b)

setsize!(y,(15,30))

@test size(y) == size(y.array)

@test size(y) == (15,30)
y .= 0
@test y == y.array

y = CachedArray(Float64,2)

@test isa(y, AbstractArray{Float64,2})

@test isa(y, CachedArray{Float64,2,Array{Float64,2}})

@test size(y) == (0,0)

setsize!(y,(15,30))

@test size(y) == (15,30)

@test isa(CachedVector(Int),CachedVector{Int,Vector{Int}})

@test isa(CachedMatrix(Int),CachedMatrix{Int,Matrix{Int}})

end # module
