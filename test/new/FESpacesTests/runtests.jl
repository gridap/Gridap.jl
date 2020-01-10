module FESpacesTests

using Test

@testset "FEFunctions" begin include("FEFunctionsTests.jl") end

@testset "FESpacesInterfaces" begin include("FESpacesInterfacesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "SingleFieldFESpaces" begin include("SingleFieldFESpacesTests.jl") end

@testset "ConformingFESpaces" begin include("ConformingFESpacesTests.jl") end

end # module
