module IntegrationTests

using Test

@testset "Quadratures" begin include("QuadraturesTests.jl") end

@testset "DuffyQuadratures" begin include("DuffyQuadraturesTests.jl") end

@testset "QuadraturesFactories" begin include("QuadratureFactoriesTests.jl") end

@testset "CellQuadratures" begin include("CellQuadraturesTests.jl") end

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "CellIntegration" begin include("CellIntegrationTests.jl") end

@testset "BoundaryTriangulations" begin include("BoundaryTriangulationsTests.jl") end

@testset "NormalVectors" begin include("NormalVectorsTests.jl") end

@testset "SkeletonTriangulations" begin include("SkeletonTriangulationsTests.jl") end

@testset "BoundaryCellFields" begin include("BoundaryCellFieldsTests.jl") end

@testset "SkeletonCellFields" begin include("SkeletonCellFieldsTests.jl") end

end # module
