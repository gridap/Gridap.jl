module FieldsTestsAll

using Test

include("FieldsMocks.jl")

include("../CellValuesTests/CellValuesMocks.jl")

include("CellFieldsMocks.jl")

@testset "FieldValues" begin include("FieldValuesTests.jl") end

@testset "FieldsMocks" begin include("FieldsMocksTests.jl") end

@testset "Fields" begin include("FieldsTests.jl") end

@testset "FieldsOperations" begin include("FieldsOperationsTests.jl") end

@testset "AnalyticalFields" begin include("AnalyticalFieldsTests.jl") end

@testset "CellFieldsMocks" begin include("CellFieldsMocksTests.jl") end

@testset "CellFieldsOperations" begin include("CellFieldsOperationsTests.jl") end

@testset "ConstantCellFields" begin include("ConstantCellFieldsTests.jl") end

end # module
