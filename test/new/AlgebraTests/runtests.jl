module AlgebraTests

using Test

@testset "LinearSolvers" begin include("LinearSolversTests.jl") end

@testset "NonLinearOperators" begin include("NonLinearOperatorsTests.jl") end

@testset "NonLinearSolvers" begin include("NonLinearSolversTests.jl") end

@testset "NLSolvers" begin include("NLSolversTests.jl") end

end # module
