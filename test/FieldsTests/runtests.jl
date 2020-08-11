module FieldsTests

using Test

@testset "MockFields" begin include("MockFieldsTests.jl") end

@testset "FieldInterface" begin include("FieldInterfaceTests.jl") end

@testset "FieldApply" begin include("FieldApplyTests.jl") end

@testset "ConstantFields" begin include("ConstantFieldsTests.jl") end

@testset "Homothecies" begin include("HomotheciesTests.jl") end

@testset "AffineMapsTests" begin include("AffineMapsTests.jl") end

@testset "FieldArrays" begin include("FieldArraysTests.jl") end

@testset "Compose" begin include("ComposeTests.jl") end

@testset "Lincomb" begin include("LincombTests.jl") end

@testset "Attachmap" begin include("AttachmapTests.jl") end

@testset "Integrate" begin include("IntegrateTests.jl") end

@testset "FieldOperations" begin include("FieldOperationsTests.jl") end

@testset "FieldOperationsNew" begin include("FieldOperationsTestsNew.jl") end

@testset "InsertInBlocks" begin include("InsertInBlocksTests.jl") end

@testset "DiffOperators" begin include("DiffOperatorsTests.jl") end

end # module
