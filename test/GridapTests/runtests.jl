module GridapTests

using Test

@testset "Poisson" begin include("PoissonTests.jl") end

@testset "PoissonDG" begin include("PoissonDGTests.jl") end

@testset "PLaplacian" begin include("PLaplacianTests.jl") end

@testset "StokesTaylorHood" begin include("StokesTaylorHoodTests.jl") end

@testset "StokesDG" begin include("StokesDGTests.jl") end

@testset "StokesNitsche" begin include("StokesNitscheTests.jl") end

@testset "Darcy" begin include("DarcyTests.jl") end

@testset "SurfaceCoupling" begin include("SurfaceCouplingTests.jl") end

@testset "IsotropicDamage" begin include("IsotropicDamageTests.jl") end

@testset "PhysicalPoisson" begin include("PhysicalPoissonTests.jl") end

end # module
