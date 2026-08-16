module DifferentialFormsTests

using Test

@testset "DifferentialFormValues" begin include("DifferentialFormValuesTests.jl") end

@testset "MetricHodge" begin include("MetricHodgeTests.jl") end

@testset "PullbackKoszul" begin include("PullbackKoszulTests.jl") end

@testset "FormBridge" begin include("FormBridgeTests.jl") end

end # module
