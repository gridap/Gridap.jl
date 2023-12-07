module TransientFEProblemsTests

using Test

@testset "HeatEquationScalar" begin include("TransientFEProblemsTests/HeatEquationScalarTests.jl") end

@testset "HeatEquationVector" begin include("TransientFEProblemsTests/HeatEquationVectorTests.jl") end

@testset "HeatEquationMultiField" begin include("TransientFEProblemsTests/HeatEquationMultiFieldTests.jl") end

@testset "HeatEquationNeumann" begin include("TransientFEProblemsTests/HeatEquationNeumannTests.jl") end

@testset "HeatEquationDG" begin include("TransientFEProblemsTests/HeatEquationDGTests.jl") end

@testset "StokesEquation" begin include("TransientFEProblemsTests/StokesEquationTests.jl") end

@testset "FreeSurfacePotentialFlow" begin include("TransientFEProblemsTests/FreeSurfacePotentialFlowTests.jl") end

@testset "SecondOrderEquation" begin include("TransientFEProblemsTests/SecondOrderEquationTests.jl") end

end # module TransientFEProblemsTests
