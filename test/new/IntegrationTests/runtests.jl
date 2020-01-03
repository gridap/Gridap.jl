module IntegrationTests

using Test

@testset "Quadratures" begin include("QuadraturesTests.jl") end

@testset "TensorProductQuadratures" begin include("TensorProductQuadraturesTests.jl") end

@testset "DuffyQuadratures" begin include("DuffyQuadraturesTests.jl") end

end # module
