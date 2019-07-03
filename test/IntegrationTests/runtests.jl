module IntegrationTests

using Test

@testset "Quadratures" begin include("QuadraturesTests.jl") end

@testset "CellQuadratures" begin include("CellQuadraturesTests.jl") end

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "CellIntegration" begin include("CellIntegrationTests.jl") end

@testset "BoundaryTriangulations" begin include("BoundaryTriangulationsTests.jl") end

@testset "BoundaryCellFields" begin include("BoundaryCellFieldsTests.jl") end

end # module
