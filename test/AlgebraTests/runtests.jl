module AlgebraTests

using Test

include("NonLinearOperatorMocks.jl")

@testset "LinearSolvers" begin include("LinearSolversTests.jl") end

@testset "NonLinearSolvers" begin include("NonLinearSolversTests.jl") end

end # module AlgebraTests
