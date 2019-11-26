module GeometryTests

using Test

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "ConformingTriangulations" begin include("ConformingTriangulationsTests.jl") end

@testset "GridGraphs" begin include("GridGraphsTests.jl") end

@testset "UnstructuredGrids" begin include("UnstructuredGridsTests.jl") end


end # module
