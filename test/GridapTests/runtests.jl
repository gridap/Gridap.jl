module GridapTests

using Test

@time @testset "Poisson" begin include("PoissonTests.jl") end

@time @testset "PoissonDG" begin include("PoissonDGTests.jl") end

# @time @testset "Poisson1D" begin include("Poisson1DTests.jl") end

# @time @testset "PLaplacian" begin include("PLaplacianTests.jl") end

# @time @testset "PLaplacianWithAutodiff" begin include("PLaplacianWithAutodiffTests.jl") end

# @time @testset "StokesTaylorHood" begin include("StokesTaylorHoodTests.jl") end

# @time @testset "StokesDG" begin include("StokesDGTests.jl") end

# @time @testset "StokesNitsche" begin include("StokesNitscheTests.jl") end

# @time @testset "Darcy" begin include("DarcyTests.jl") end

# @time @testset "PeriodicDarcy" begin include("PeriodicDarcyTests.jl") end

# @time @testset "SurfaceCoupling" begin include("SurfaceCouplingTests.jl") end

# @time @testset "IsotropicDamage" begin include("IsotropicDamageTests.jl") end

# @time @testset "PhysicalPoisson" begin include("PhysicalPoissonTests.jl") end

end # module
