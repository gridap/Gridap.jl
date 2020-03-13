module FESpacesTests

using Test

@testset "CellBases" begin include("CellBasesTests.jl") end

@testset "Law" begin include("LawTests.jl") end

@testset "FEFunctions" begin include("FEFunctionsTests.jl") end

@testset "FESpacesInterfaces" begin include("FESpacesInterfacesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "FEOperators" begin include("FEOperatorsTests.jl") end

@testset "SingleFieldFESpaces" begin include("SingleFieldFESpacesTests.jl") end

@testset "SingleFieldFEFunctions" begin include("SingleFieldFEFunctionsTests.jl") end

@testset "TrialFESpaces" begin include("TrialFESpacesTests.jl") end

@testset "SparseMatrixAssemblers" begin include("SparseMatrixAssemblersTests.jl") end

@testset "UnconstrainedFESpaces" begin include("UnconstrainedFESpacesTests.jl") end

@testset "ConformingFESpaces" begin include("ConformingFESpacesTests.jl") end

@testset "DivConformingFESpaces" begin include("DivConformingFESpacesTests.jl") end

@testset "CurlConformingFESpaces" begin include("CurlConformingFESpacesTests.jl") end

@testset "DiscontinuousFESpaces" begin include("DiscontinuousFESpacesTests.jl") end

@testset "FETerms" begin include("FETermsTests.jl") end

@testset "CellKernels" begin include("CellKernelsTests.jl") end

@testset "AffineFEOperators" begin include("AffineFEOperatorsTests.jl") end

@testset "FEOperatorsFromTerms" begin include("FEOperatorsFromTermsTests.jl") end

@testset "FESolvers" begin include("FESolversTests.jl") end

@testset "FESpacesWithLastDofRemoved" begin include("FESpacesWithLastDofRemovedTests.jl") end

@testset "ZeroMeanFESpaces" begin include("ZeroMeanFESpacesTests.jl") end

@testset "CLagrangianFESpaces" begin include("CLagrangianFESpacesTests.jl") end

@testset "DirichletFESpaces" begin include("DirichletFESpacesTests.jl") end

@testset "ExtendedFESpaces" begin include("ExtendedFESpacesTests.jl") end

@testset "FESpaceFactories" begin include("FESpaceFactoriesTests.jl") end

@testset "StateLaws" begin include("StateLawsTests.jl") end

end # module
