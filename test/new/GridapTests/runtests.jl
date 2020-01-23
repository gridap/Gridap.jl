module GridapTests

using Test

@testset "Poisson" begin include("PoissonTests.jl") end

@testset "PoissonDG" begin include("PoissonDGTests.jl") end

end # module
