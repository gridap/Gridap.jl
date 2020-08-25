module ArraysTests

using Test

@testset "Interfaces" begin include("InterfaceTests.jl") end

@testset "BlockArraysCoo" begin include("BlockArraysCooTests.jl") end

@testset "VectorsOfBlockArrayCoo" begin include("VectorsOfBlockArrayCooTests.jl") end

@testset "CachedArrays" begin include("CachedArraysTests.jl") end

@testset "Kernels" begin include("KernelsTests.jl") end

@testset "Apply" begin include("ApplyTests.jl") end

@testset "CompressedArrays" begin include("CompressedArraysTests.jl") end

@testset "LocalToGlobalArrays" begin include("LocalToGlobalArraysTests.jl") end

@testset "LocalToGlobalPosNegArrays" begin include("LocalToGlobalPosNegArraysTests.jl") end

@testset "FilteredArraysTests" begin include("FilteredArraysTests.jl") end

@testset "Tables" begin include("TablesTests.jl") end

@testset "Reindex" begin include("ReindexTests.jl") end

@testset "IdentityVectors" begin include("IdentityVectorsTests.jl") end

@testset "SubVectors" begin include("SubVectorsTests.jl") end

@testset "ArrayPairs" begin include("ArrayPairsTests.jl") end

@testset "AppendedArrays" begin include("AppendedArraysTests.jl") end

@testset "AutodiffTests" begin include("AutodiffTests.jl") end

end # module
