module GridapRunTests

using Test

@time @testset "Helpers" begin include("HelpersTests/runtests.jl") end

@time @testset "Io" begin include("IoTests/runtests.jl") end

@time @testset "Algebra" begin include("AlgebraTests/runtests.jl") end

@time @testset "Arrays" begin include("ArraysTests/runtests.jl") end

@time @testset "TensorValues" begin include("TensorValuesTests/runtests.jl") end

@time @testset "Fields" begin include("FieldsTests/runtests.jl") end

@time @testset "Polynomials" begin include("PolynomialsTests/runtests.jl") end

@time @testset "Integration" begin include("IntegrationTests/runtests.jl") end

@time @testset "ReferenceFEs" begin include("ReferenceFEsTests/runtests.jl") end

@time @testset "Geometry" begin include("GeometryTests/runtests.jl") end

@time @testset "CellData" begin include("CellDataTests/runtests.jl") end

@time @testset "Visualization" begin include("VisualizationTests/runtests.jl") end

@time @testset "FESpaces (1/2)" begin include("FESpacesTests/runtests_1.jl") end

@time @testset "FESpaces (2/2)" begin include("FESpacesTests/runtests_2.jl") end

@time @testset "MultiField" begin include("MultiFieldTests/runtests.jl") end

end # module
