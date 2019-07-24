module GeometryTests

using Test

@testset "GridGraphs" begin include("GridGraphsTests.jl") end

@testset "Grids" begin include("GridsTests.jl") end

@testset "FaceLabels" begin include("FaceLabelsTests.jl") end

@testset "DiscreteModels" begin include("DiscreteModelsTests.jl") end

@testset "UnstructuredGeometry" begin include("UnstructuredGeometryTests.jl") end

@testset "GeometryWrappers" begin include("GeometryWrappersTests.jl") end

@testset "CartesianGeometry" begin include("CartesianGeometryTests.jl") end

@testset "GridPortions" begin include("GridPortionsTests.jl") end

@testset "BoundaryGrids" begin include("BoundaryGridsTests.jl") end

@testset "SkeletonGrids" begin include("SkeletonGridsTests.jl") end

end # module
