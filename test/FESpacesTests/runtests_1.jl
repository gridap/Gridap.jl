module FESpacesTests1

using Test

@testset "ConformingFESpaces" begin include("ConformingFESpacesTests.jl") end

@testset "FESpacesInterfaces" begin include("FESpacesInterfacesTests.jl") end

@testset "SingleFieldFESpaces" begin include("SingleFieldFESpacesTests.jl") end

@testset "TrialFESpaces" begin include("TrialFESpacesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "SparseMatrixAssemblers" begin include("SparseMatrixAssemblersTests.jl") end

@testset "FEOperators" begin include("FEOperatorsTests.jl") end

@testset "AffineFEOperators" begin include("AffineFEOperatorsTests.jl") end

@testset "FEOperatorsFromWeakForm" begin include("FEOperatorsFromWeakFormTests.jl") end

@testset "FESolvers" begin include("FESolversTests.jl") end

@testset "DiscontinuousFESpaces" begin include("DiscontinuousFESpacesTests.jl") end

@testset "DivConformingFESpaces" begin include("DivConformingFESpacesTests.jl") end

@testset "CurlConformingFESpaces" begin include("CurlConformingFESpacesTests.jl") end

end # module
