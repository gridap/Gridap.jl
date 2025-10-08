module FESpacesTests2

using Test

@testset "FESpacesWithConstantFixed" begin include("FESpacesWithConstantFixedTests.jl") end

@testset "ZeroMeanFESpaces" begin include("ZeroMeanFESpacesTests.jl") end

@testset "CLagrangianFESpaces" begin include("CLagrangianFESpacesTests.jl") end

@testset "DirichletFESpaces" begin include("DirichletFESpacesTests.jl") end

@testset "ExtendedFESpaces" begin include("ExtendedFESpacesTests.jl") end

@testset "FESpacesWithLinearConstraints" begin include("FESpacesWithLinearConstraintsTests.jl") end

@testset "FEAutodiff" begin include("FEAutodiffTests.jl") end

@testset "CDLagrangianFESpaces" begin include("CDLagrangianFESpacesTests.jl") end

@testset "PhysicalFESpaces" begin include("PhysicalFESpacesTests.jl") end

@testset "FESpaceFactories" begin include("FESpaceFactoriesTests.jl") end

@testset "AppendedTriangulations" begin include("AppendedTriangulationsTests.jl") end

@testset "DiscreteModelWithFEMapsTests" begin include("DiscreteModelWithFEMapsTests.jl") end

@testset "ConstantFESpaceTests" begin include("ConstantFESpaceTests.jl") end

@testset "PatchAssemblersTests" begin include("PatchAssemblersTests.jl") end

@testset "PolytopalFESpacesTests" begin include("PolytopalFESpacesTests.jl") end

end # module
