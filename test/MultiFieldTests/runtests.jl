module MultiFieldTests

using Test

@testset "MultiCellArrays" begin include("MultiCellArraysTests.jl") end
@testset "MultiCellMaps" begin include("MultiCellMapsTests.jl") end
@testset "MultiAssemblers" begin include("MultiAssemblersTests.jl") end
@testset "MultiFESpaces" begin include("MultiFESpacesTests.jl") end
@testset "MultiFEFunctions" begin include("MultiFEFunctionsTests.jl") end
@testset "MultiFEOperators" begin include("MultiFEOperatorsTests.jl") end

end # module MultiFieldTests
