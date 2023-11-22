module ODESolversAllTests

using Test

using Gridap
using Gridap.ODEs

@time @testset "ThetaMethodSolversFamilyTests" begin include("ThetaMethodSolversFamilyTests.jl") end

@time @testset "TableausTests" begin include("TableausTests.jl") end

end # module ODESolversAllTests
