module ArraysTests

using Test

@testset "Interfaces" begin include("InterfaceTests.jl") end

@testset "CachedArrays" begin include("CachedArraysTests.jl") end

@testset "Kernels" begin include("KernelsTests.jl") end

@testset "Apply" begin include("ApplyTests.jl") end

@testset "CompressedArrays" begin include("CompressedArraysTests.jl") end

@testset "LocalToGlobalArrays" begin include("LocalToGlobalArraysTests.jl") end

@testset "Tables" begin include("TablesTests.jl") end

@testset "Reindex" begin include("ReindexTests.jl") end

end # module
