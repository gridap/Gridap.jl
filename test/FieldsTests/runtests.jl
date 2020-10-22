module FieldsTests

using Test

@testset "Fields" begin include("FieldInterfacesTests.jl") end

#@testset "FieldOperations" begin include("FieldOperationsTests.jl") end

@testset "FieldArrays" begin include("FieldArraysTests.jl") end

@testset "FieldArraysOperations" begin include("FieldArraysOperationsTests.jl") end

#@testset "LazyFieldArrays" begin include("LazyFieldArraysTests.jl") end

end
