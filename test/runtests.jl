using Numa
using Test

@testset "Numa.jl" begin
  @time @testset "CellArrays" begin include("cellarraystests.jl") end
  @time @testset "CellQuadratures" begin include("cellquadraturestests.jl") end
  @time @testset "Polynomials" begin include("polynomialtests.jl") end
  @time @testset "Polytopes" begin include("polytopetests.jl") end
  @time @testset "RefFEs" begin include("reffetests.jl") end
  @time @testset "Meshes" begin include("meshtests.jl") end
  @time @testset "FESpaces" begin include("fespacetests.jl") end
  @time @testset "BilinearForms" begin include("bilinearformtests.jl") end
end
