module GridapTests

using Test

@testset "Poisson" begin include("PoissonTests.jl") end

@testset "PoissonDG" begin include("PoissonDGTests.jl") end

@testset "PLaplacian" begin include("PLaplacianTests.jl") end

@testset "StokesTaylorHood" begin include("StokesTaylorHoodTests.jl") end

@testset "StokesDG" begin include("StokesDGTests.jl") end

@testset "Darcy" begin include("DarcyTests.jl") end

end # module
