module MultiFieldTests

using Test

@testset "MultiCellArrays" begin include("MultiCellArraysTests.jl") end

@testset "MultiFieldCellBases" begin include("MultiFieldCellBasesTests.jl") end

end # module
