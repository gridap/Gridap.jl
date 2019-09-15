module FESpacesTestsAll

using Test

@testset "FESpaces" begin include("FESpacesTests.jl") end

@testset "ConformingFESpaces" begin include("ConformingFESpacesTests.jl") end

@testset "CLagrangianFESpaces" begin include("CLagrangianFESpacesTests.jl") end

@testset "DLagrangianFESpaces" begin include("DLagrangianFESpacesTests.jl") end

@testset "DiscFESpaces" begin include("DiscFESpacesTests.jl") end

@testset "RaviartThomasFESpacesTests" begin include("RaviartThomasFESpacesTests.jl") end

@testset "ConstrainedFESpaces" begin include("ConstrainedFESpacesTests.jl") end

@testset "FEFunctions" begin include("FEFunctionsTests.jl") end

@testset "FEBases" begin include("FEBasesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "FEterms" begin include("FETermsTests.jl") end

@testset "FEOperators" begin include("FEOperatorsTests.jl") end

@testset "DGFEOperators" begin include("DGFEOperatorsTests.jl") end

@testset "VectorValuedFEOperators" begin include("VectorValuedFEOperatorsTests.jl") end

@testset "NonLinearFEOperators" begin include("NonLinearFEOperatorsTests.jl") end

@testset "VectorValuedNonLinearFEOperators" begin include("VectorValuedNonLinearFEOperatorsTests.jl") end

end # module
