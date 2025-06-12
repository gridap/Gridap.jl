module FieldsTests

using Test

@testset "Fields" begin include("FieldInterfacesTests.jl") end

@testset "FieldOperations" begin include("FieldOperationsTests.jl") end

@testset "FieldArrays" begin include("FieldArraysTests.jl") end

@testset "DiffOperators" begin include("DiffOperatorsTests.jl") end

@testset "AffineMaps" begin include("AffineMapsTests.jl") end

@testset "FieldArraysOperations" begin include("FieldArraysOperationsTests.jl") end

@testset "LazyFieldArrays" begin include("LazyFieldArraysTests.jl") end

@testset "FieldArrayBlocks" begin include("FieldArrayBlocksTests.jl") end

@testset "InverseFields" begin include("InverseFieldsTests.jl") end

@testset "DensifyInnerMostBlockLevelMapsTests" begin include("DensifyInnerMostBlockLevelMapsTests.jl") end

end
