module NewFieldsTests

using Test

@testset "Fields" begin include("FieldInterfacesTests.jl") end

# @testset "MockFields" begin include("MockFieldsTests.jl") end

# @testset "ConstantFields" begin include("ConstantFieldsTests.jl") end

# @testset "FunctionFields" begin include("FunctionFieldsTests.jl") end

# @testset "LinearCombinationFields" begin include("LinearCombinationFieldsTests.jl") end

@testset "FieldOperations" begin include("FieldOperationsTests.jl") end

# @testset "FieldArrays" begin include("FieldArraysTests.jl") end

end
