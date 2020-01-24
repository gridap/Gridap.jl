module InferenceTests

using Test
using Gridap.Inference
import Gridap.Inference: testargs

bar(x,y) = x*y
foo(x) = sqrt(x-1)
hola(x,y) = sqrt(x-1)*y

testargs(::typeof(foo),T::DataType) = (one(T),)
testargs(::typeof(hola),Tx,Ty) = (one(Tx),testvalue(Ty))

@test 0 == testvalue(Int)
@test (0,0.0) == testvalues(Int,Float64)

@test return_type(sqrt, Int) == Float64
@test return_type(foo, Int) == Float64
@test return_type(bar, Int, Int) == Int
@test return_type(bar, Float64, Int) == Float64
@test return_type(bar, Matrix{Float64}, Int) == Matrix{Float64}
@test return_type(hola, Int, Vector{Int}) == Vector{Float64}

@test return_type_broadcast(bar, Matrix{Float64}, Int) == Matrix{Float64}
@test return_type_broadcast(foo, Matrix{Float64}) == Matrix{Float64}
@test return_type_broadcast(foo, Matrix{Float64}) == Matrix{Float64}
@test return_type_broadcast(bar, Float64, Int) == Float64
@test return_type_broadcast(foo, Float64) == Float64
@test return_type_broadcast(bar, Int, Int) == Int
@test return_type_broadcast(bar, Int, Vector{Int}) == Vector{Int}
@test return_type_broadcast(hola, Int, Vector{Int}) == Vector{Float64}
@test return_type_broadcast(hola, Int, Int) == Float64

end # module
