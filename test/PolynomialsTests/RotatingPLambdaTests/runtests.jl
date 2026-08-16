module RotatingPLambdaTests

using Test

@testset "PΛRotations" begin include("PΛRotationsTests.jl") end

@testset "PΛTrimmedRotations" begin include("PΛTrimmedRotationsTests.jl") end

@testset "PΛRotationPullback" begin include("PΛRotationPullbackTests.jl") end

@testset "PΛBasisGradients" begin include("PΛBasisGradientTests.jl") end

@testset "RotatingPΛBases" begin include("RotatingPΛBasesTests.jl") end

@testset "TrimmedPΛBases" begin include("TrimmedPΛBasesTests.jl") end

end # module
