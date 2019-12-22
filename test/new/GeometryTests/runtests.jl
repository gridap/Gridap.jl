module GeometryTests

using Test

@testset "GridTopologies" begin include("GridTopologiesTests.jl") end

@testset "UnstructuredGridTopologies" begin include("UnstructuredGridTopologiesTests.jl") end

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "Grids" begin include("GridsTests.jl") end

@testset "DiscreteModels" begin include("DiscreteModelsTests.jl") end

@testset "FaceLabelings" begin include("FaceLabelingsTests.jl") end

@testset "UnstructuredGrids" begin include("UnstructuredGridsTests.jl") end

@testset "CartesianGrids" begin include("CartesianGridsTests.jl") end

@testset "UnstructuredDiscreteModels" begin include("UnstructuredDiscreteModelsTests.jl") end

@testset "CartesianDiscreteModels" begin include("CartesianDiscreteModelsTests.jl") end

end # module
