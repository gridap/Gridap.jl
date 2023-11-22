module ODESolversAllTests

using Test

using Gridap
using Gridap.ODEs

@time @testset "ThetaMethodSolversFamilyTests.jl" begin include("ThetaMethodSolversFamilyTests.jl") end

end # module ODESolversAllTests
