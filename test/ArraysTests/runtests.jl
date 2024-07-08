module ArraysTests

using Test

@testset "Interfaces" begin include("InterfaceTests.jl") end

@testset "CachedArrays" begin include("CachedArraysTests.jl") end

@testset "Maps" begin include("MapsTests.jl") end

@testset "AlgebraMaps" begin include("AlgebraMapsTests.jl") end

@testset "LazyArrays" begin include("LazyArraysTests.jl") end

@testset "CompressedArrays" begin include("CompressedArraysTests.jl") end

@testset "FilteredArraysTests" begin include("FilteredArraysTests.jl") end

@testset "Tables" begin include("TablesTests.jl") end

@testset "Reindex" begin include("ReindexTests.jl") end

@testset "KeyToValMapsTests" begin include("KeyToValMapsTests.jl") end

@testset "PosNegReindex" begin include("PosNegReindexTests.jl") end

@testset "IdentityVectors" begin include("IdentityVectorsTests.jl") end

@testset "ArrayPairs" begin include("ArrayPairsTests.jl") end

@testset "AppendedArrays" begin include("AppendedArraysTests.jl") end

@testset "AutodiffTests" begin include("AutodiffTests.jl") end

@testset "VectorWithEntryRemovedTests" begin include("VectorWithEntryRemovedTests.jl") end

@testset "VectorWithEntryInsertedTests" begin include("VectorWithEntryInsertedTests.jl") end

@testset "PrintOpTrees" begin include("PrintOpTreesTests.jl") end

end # module
