module TransientFEProblemsTests

using Test

@testset "Order1FETests" begin include("TransientFEProblemsTests/Order1FETests.jl") end

@testset "Order2FETests" begin include("TransientFEProblemsTests/Order2FETests.jl") end

@testset "HeatEquationVectorTests" begin include("TransientFEProblemsTests/HeatEquationVectorTests.jl") end

@testset "HeatEquationMultiFieldTests" begin include("TransientFEProblemsTests/HeatEquationMultiFieldTests.jl") end

# @testset "HeatEquationBoundaryTests" begin include("TransientFEProblemsTests/HeatEquationBoundaryTests.jl") end

# @testset "HeatEquationDGTests" begin include("TransientFEProblemsTests/HeatEquationDGTests.jl") end

# @testset "StokesEquationTests" begin include("TransientFEProblemsTests/StokesEquationTests.jl") end

# @testset "StokesEquationAutoDiffTests" begin include("TransientFEProblemsTests/StokesEquationAutoDiffTests.jl") end

# @testset "FreeSurfacePotentialFlowTests" begin include("TransientFEProblemsTests/FreeSurfacePotentialFlowTests.jl") end

end # module TransientFEProblemsTests
