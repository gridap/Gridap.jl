module FieldsTests

using Test

include("FieldsMocks.jl")

@testset "FieldValues" begin include("FieldValuesTests.jl") end

@testset "FieldsMocks" begin include("FieldsMocksTests.jl") end

@testset "Fields" begin include("FieldsTests.jl") end

@testset "FieldsOperations" begin include("FieldsOperationsTests.jl") end

end # module
