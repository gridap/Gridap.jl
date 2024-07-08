module GeometryTests

using Test

@testset "Grids" begin include("GridsTests.jl") end

@testset "UnstructuredGrids" begin include("UnstructuredGridsTests.jl") end

@testset "CartesianGrids" begin include("CartesianGridsTests.jl") end

@testset "GridTopologies" begin include("GridTopologiesTests.jl") end

@testset "UnstructuredGridTopologies" begin include("UnstructuredGridTopologiesTests.jl") end

@testset "FaceLabelings" begin include("FaceLabelingsTests.jl") end

@testset "DiscreteModels" begin include("DiscreteModelsTests.jl") end

@testset "MappedDiscreteModels" begin include("MappedDiscreteModelsTests.jl") end

@testset "UnstructuredDiscreteModels" begin include("UnstructuredDiscreteModelsTests.jl") end

@testset "CartesianDiscreteModels" begin include("CartesianDiscreteModelsTests.jl") end

@testset "PeriodicBC" begin include("PeriodicBCTests.jl") end

@testset "GridPortions" begin include("GridPortionsTests.jl") end

@testset "DiscreteModelPortions" begin include("DiscreteModelPortionsTests.jl") end

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "BoundaryTriangulations" begin include("BoundaryTriangulationsTests.jl") end

@testset "SkeletonTriangulations" begin include("SkeletonTriangulationsTests.jl") end

@testset "AppendedTriangulations" begin include("AppendedTriangulationsTests.jl") end

@testset "CompressedCellArrays" begin include("CompressedCellArraysTests.jl") end

end # module
