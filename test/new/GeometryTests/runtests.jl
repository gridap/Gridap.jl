module GeometryTests

using Test

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "ConformingTriangulations" begin include("ConformingTriangulationsTests.jl") end

@testset "DiscreteModels" begin include("DiscreteModelsTests.jl") end

@testset "UnstructuredGrids" begin include("UnstructuredGridsTests.jl") end

@testset "CartesianGrids" begin include("CartesianGridsTests.jl") end

@testset "UnstructuredDiscreteModels" begin include("UnstructuredDiscreteModelsTests.jl") end

end # module
