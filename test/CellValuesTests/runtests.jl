module CellValuesTests

using Test

include("CellValuesMocks.jl")

include("MapsMocks.jl")

@testset "CellValuesMocks" begin include("CellValuesMocksTests.jl") end

@testset "MapsMocks" begin include("MapsMocksTests.jl") end

@testset "CellNumbers" begin include("CellNumbersTests.jl") end

@testset "CellArrays" begin include("CellArraysTests.jl") end

@testset "CellMaps" begin include("CellMapsTests.jl") end

@testset "Kernels" begin include("KernelsTests.jl") end

@testset "CellNumberApply" begin include("CellNumberApplyTests.jl") end

@testset "CellArrayApply" begin include("CellArrayApplyTests.jl") end

@testset "MapApply" begin include("MapApplyTests.jl") end

@testset "CellMapApply" begin include("CellMapApplyTests.jl") end

@testset "CellValuesOperations" begin include("CellValuesOperationsTests.jl") end

@testset "CellValuesGallery" begin include("CellValuesGalleryTests.jl") end

@testset "CellValuesAppend" begin include("CellValuesAppendTests.jl") end

@testset "CellValuesReindex" begin include("CellValuesReindexTests.jl") end

@testset "FindLocalIndex" begin include("FindLocalIndexTests.jl") end

@testset "ConstantCellValues" begin include("ConstantCellValuesTests.jl") end


end # module

