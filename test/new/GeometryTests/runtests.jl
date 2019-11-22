module GeometryTests

using Test

@testset "Triangulations" begin include("TriangulationsTests.jl") end

@testset "ConformingTriangulations" begin include("ConformingTriangulationsTests.jl") end

end # module
