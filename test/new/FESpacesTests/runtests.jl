module FESpacesTests

using Test

@testset "CellBases" begin include("CellBasesTests.jl") end

@testset "FEFunctions" begin include("FEFunctionsTests.jl") end

@testset "FESpacesInterfaces" begin include("FESpacesInterfacesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "SingleFieldFESpaces" begin include("SingleFieldFESpacesTests.jl") end

@testset "SingleFieldFEFunctions" begin include("SingleFieldFEFunctionsTests.jl") end

@testset "TrialFESpaces" begin include("TrialFESpacesTests.jl") end

@testset "SingleFieldAssemblers" begin include("SingleFieldAssemblersTests.jl") end

@testset "UnconstrainedFESpaces" begin include("UnconstrainedFESpacesTests.jl") end

@testset "ConformingFESpaces" begin include("ConformingFESpacesTests.jl") end

end # module
