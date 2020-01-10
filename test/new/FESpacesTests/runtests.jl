module FESpacesTests

include("../../../src/new/FESpaces/FESpaces.jl")

using Test

@testset "FEFunctions" begin include("FEFunctionsTests.jl") end

@testset "FESpacesInterfaces" begin include("FESpacesInterfacesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "SingleFieldFESpaces" begin include("SingleFieldFESpacesTests.jl") end

end # module
