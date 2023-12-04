module ODESolversAllTests

using Test

@time @testset "TableausTests" begin include("TableausTests.jl") end

@time @testset "Order1Tests" begin include("Order1Tests.jl") end

@time @testset "Order2Tests" begin include("Order2Tests.jl") end

end # module ODESolversAllTests
