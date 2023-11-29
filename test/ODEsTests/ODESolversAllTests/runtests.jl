module ODESolversAllTests

using Test

using Gridap
using Gridap.ODEs

# @time @testset "ThetaMethodFamilyTests" begin include("ThetaMethodFamilyTests.jl") end

# @time @testset "TableausTests" begin include("TableausTests.jl") end

@time @testset "RungeKuttaTests" begin include("RungeKuttaTests.jl") end

end # module ODESolversAllTests
