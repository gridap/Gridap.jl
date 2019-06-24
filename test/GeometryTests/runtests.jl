module GeometryTestsAll

using Test

@testset "Geometry" begin include("GeometryTests.jl") end

@testset "GridPortions" begin include("GridPortionsTests.jl") end

end # module
