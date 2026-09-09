module ConvergenceTests

# IMPORTANT: These tests are NOT run by default. They are quite expensive and redundant. 

using Test

@testset "Hybrid Methods" begin
  @time @testset "Poisson - HDG" begin include("HDG.jl") end
  @time @testset "Poisson - HHO" begin include("HHO.jl") end
  @time @testset "Poisson - HHO (mixed order)" begin include("HHOMixed.jl") end
  @time @testset "Poisson - HDG polytopal" begin include("HDGPolytopal.jl") end
  @time @testset "Poisson - HHO polytopal" begin include("HHOPolytopal.jl") end
  @time @testset "Poisson - HHO polytopal (mixed order)" begin include("HHOMixedPolytopal.jl") end
  @time @testset "Elasticity - HHO (mixed order)" begin include("HHOMixedElasticity.jl") end
end

@testset "Fourth-order plate elements" begin
  @time @testset "Biharmonic - Argyris" begin include("Argyris.jl") end
  @time @testset "Biharmonic - Morley" begin include("Morley.jl") end
  @time @testset "Kirchhoff plate - HHJ" begin include("HHJ.jl") end
end

@testset "Mixed elasticity" begin
  @time @testset "Hellinger-Reissner - Arnold-Winther" begin include("ArnoldWinther.jl") end
end

@testset "Darcy-Stokes" begin
  @time @testset "Mardal-Tai-Winther, 2D and 3D" begin include("MardalTaiWinther.jl") end
end

end