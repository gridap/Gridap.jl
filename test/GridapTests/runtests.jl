module GridapTests

using Test

@testset "Conforming Methods" begin
  @time @testset "PLaplacian" begin include("PLaplacianTests.jl") end
  @time @testset "Poisson1D" begin include("Poisson1DTests.jl") end
  @time @testset "PoissonLagrangeMultiplier" begin include("PoissonLagrangeMultiplierTests.jl") end
  @time @testset "Poisson" begin include("PoissonTests.jl") end
  @time @testset "Biharmonic" begin include("BiharmonicTests.jl") end
  @time @testset "Darcy" begin include("DarcyTests.jl") end
  @time @testset "PeriodicDarcy" begin include("PeriodicDarcyTests.jl") end
  @time @testset "HCurl" begin include("HCurlTests.jl") end
  @time @testset "StokesTaylorHood" begin include("StokesTaylorHoodTests.jl") end
  @time @testset "StokesNitsche" begin include("StokesNitscheTests.jl") end
  @time @testset "StokesMini" begin include("StokesMiniTests.jl") end
end

@testset "DG Methods" begin
  @time @testset "Poisson - DG" begin include("PoissonDGTests.jl") end
  @time @testset "Poisson - DG polytopal" begin include("PoissonDGPolytopalTests.jl") end
  @time @testset "Stokes - DG" begin include("StokesDGTests.jl") end
end

@testset "Hybrid Methods" begin
  @time @testset "Poisson - HDG" begin include("HDGTests.jl") end
  @time @testset "Poisson - HHO" begin include("HHOTests.jl") end
  @time @testset "Poisson - HHO (mixed order)" begin include("HHOMixedTests.jl") end
  @time @testset "Poisson - HDG polytopal" begin include("HDGPolytopalTests.jl") end
  @time @testset "Poisson - HHO polytopal" begin include("HHOPolytopalTests.jl") end
  @time @testset "Poisson - HHO polytopal (mixed order)" begin include("HHOMixedPolytopalTests.jl") end
  @time @testset "Elasticity - HHO (mixed order)" begin include("HHOMixedElasticity.jl") end
  @time @testset "Stokes - HHO (mixed order)" begin include("HHOMixedStokes.jl") end
  @time @testset "Stokes - HHO polytopal (mixed order)" begin include("HHOMixedStokesPolytopal.jl") end
end

@time @testset "SurfaceCoupling" begin include("SurfaceCouplingTests.jl") end

@time @testset "IsotropicDamage" begin include("IsotropicDamageTests.jl") end

@testset "Issue reproducers" begin
  @time @testset "Issue614" begin include("issues/issue_614.jl") end
  @time @testset "Issue689" begin include("issues/issue_689.jl") end
  @time @testset "Issue722" begin include("issues/issue_722.jl") end
  @time @testset "Issue743" begin include("issues/issue_743.jl") end
  @time @testset "Issue760" begin include("issues/issue_760.jl") end
  @time @testset "Issue770" begin include("issues/issue_770.jl") end
  @time @testset "Issue778" begin include("issues/issue_778.jl") end
  @time @testset "Issue869" begin include("issues/issue_869.jl") end
end

@time @testset "EmptyDomains" begin include("EmptyDomainsTests.jl") end

end # module
