module GridapTests

using Test

@time @testset "Poisson" begin include("PoissonTests.jl") end

@time @testset "PoissonDG" begin include("PoissonDGTests.jl") end

@time @testset "Poisson1D" begin include("Poisson1DTests.jl") end

@time @testset "PLaplacian" begin include("PLaplacianTests.jl") end

@time @testset "StokesTaylorHood" begin include("StokesTaylorHoodTests.jl") end

@time @testset "StokesDG" begin include("StokesDGTests.jl") end

@time @testset "StokesNitsche" begin include("StokesNitscheTests.jl") end

@time @testset "Darcy" begin include("DarcyTests.jl") end

@time @testset "PeriodicDarcy" begin include("PeriodicDarcyTests.jl") end

@time @testset "SurfaceCoupling" begin include("SurfaceCouplingTests.jl") end

@time @testset "IsotropicDamage" begin include("IsotropicDamageTests.jl") end

@time @testset "Biharmonic" begin include("BiharmonicTests.jl") end

@time @testset "PoissonLagrangeMultiplier" begin include("PoissonLagrangeMultiplierTests.jl") end

@time @testset "Issue614" begin include("issue_614.jl") end

end # module
