module ODESolversAllTests

using Test

@time @testset "Tableaus" begin include("ODEProblemsTests/TableausTests.jl") end

@time @testset "Order1ODE" begin include("ODEProblemsTests/Order1ODETests.jl") end

@time @testset "Order2ODE" begin include("ODEProblemsTests/Order2ODETests.jl") end

end # module ODESolversAllTests
