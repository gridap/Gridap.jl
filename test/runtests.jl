module GridapTests

using Test

@time @testset "Utils" begin include("UtilsTests/runtests.jl") end

@time @testset "CellValues" begin include("CellValuesTests/runtests.jl") end

@time @testset "Fields" begin include("FieldsTests/runtests.jl") end

@time @testset "RefFEs" begin include("RefFEsTests/runtests.jl") end

@time @testset "Integration" begin include("IntegrationTests/runtests.jl") end

@time @testset "Geometry" begin include("GeometryTests/runtests.jl") end

@time @testset "Algebra" begin include("AlgebraTests/runtests.jl") end

@time @testset "FESpaces" begin include("FESpacesTests/runtests.jl") end

@time @testset "MultiField" begin include("MultiFieldTests/runtests.jl") end

@time @testset "Visualization" begin include("VisualizationTests/runtests.jl") end

end # module
