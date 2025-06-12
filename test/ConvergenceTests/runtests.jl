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

end