module MultiFieldTests

using Test

@testset "MultiFieldCellFields" begin include("MultiFieldCellFieldsTests.jl") end

@testset "MultiFieldFESpaces" begin include("MultiFieldFESpacesTests.jl") end

@testset "MultiFieldFEFunctions" begin include("MultiFieldFEFunctionsTests.jl") end

@testset "MultiFieldSparseMatrixAssemblers" begin include("MultiFieldSparseMatrixAssemblersTests.jl") end

@testset "MultiFieldFEOperators" begin include("MultiFieldFEOperatorsTests.jl") end

@testset "MultiFieldFESpacesWithLinearConstraints" begin include("MultiFieldFESpacesWithLinearConstraintsTests.jl") end

end # module
