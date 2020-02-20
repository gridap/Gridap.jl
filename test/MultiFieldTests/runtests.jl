module MultiFieldTests

using Test

@testset "MultiFieldArrays" begin include("MultiFieldArraysTests.jl") end

@testset "MultiCellArrays" begin include("MultiCellArraysTests.jl") end

@testset "MultiFieldCellBases" begin include("MultiFieldCellBasesTests.jl") end

@testset "MultiFieldFESpaces" begin include("MultiFieldFESpacesTests.jl") end

@testset "MultiFieldFEFunctions" begin include("MultiFieldFEFunctionsTests.jl") end

@testset "MultiFieldSparseMatrixAssemblers" begin include("MultiFieldSparseMatrixAssemblersTests.jl") end

@testset "MultiFieldFEOperators" begin include("MultiFieldFEOperatorsTests.jl") end

end # module
