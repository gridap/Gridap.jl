module FieldsTests

using Test

@testset "FieldValues" begin include("FieldValuesTests.jl") end

@testset "Fields" begin include("FieldsTests.jl") end

end # module
