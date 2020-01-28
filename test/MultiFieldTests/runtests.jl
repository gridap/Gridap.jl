module MultiFieldTests

using Test

@testset "MultiCellArrays" begin include("MultiCellArraysTests.jl") end

@testset "MultiFieldCellBases" begin include("MultiFieldCellBasesTests.jl") end

@testset "MultiFieldFESpaces" begin include("MultiFieldFESpacesTests.jl") end

@testset "MultiFieldFEFunctions" begin include("MultiFieldFEFunctionsTests.jl") end

end # module
