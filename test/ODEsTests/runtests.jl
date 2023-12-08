module ODEsTests

using Test

@time @testset "TimeDerivatives" begin include("TimeDerivativesTests.jl") end

@time @testset "ODEOperators" begin include("ODEOperatorsTests.jl") end

@time @testset "ODESolvers" begin include("ODESolversTests.jl") end

@time @testset "ODEProblems" begin include("ODEProblemsTests.jl") end

@time @testset "ODESolutions" begin include("ODESolutionsTests.jl") end

@time @testset "TransientFESpaces" begin include("TransientFESpacesTests.jl") end

@time @testset "TransientCellFields" begin include("TransientCellFieldsTests.jl") end

@time @testset "TransientFEOperatorsSolutions" begin include("TransientFEOperatorsSolutionsTests.jl") end

@time @testset "TransientFEProblems" begin include("TransientFEProblemsTests.jl") end

# @time @testset "DiffEqsWrappers" begin include("_DiffEqsWrappersTests.jl") end

# include("../bench/runbenchs.jl")

end # module ODEsTests
