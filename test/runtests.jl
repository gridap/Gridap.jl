using Numa
using Test

@testset "Numa.jl" begin
  @time @testset "Quadratures" begin include("quadraturestests.jl") end
  @time @testset "CellArrays" begin include("cellarraystests.jl") end
  @time @testset "CellQuadratures" begin include("cellquadraturestests.jl") end
  @time @testset "Polynomials" begin include("polynomialtests.jl") end
  @time @testset "CellBases" begin include("cellbasestests.jl") end
  @time @testset "CellFields" begin include("cellfieldstests.jl") end
  @time @testset "CellScalarsVectorsAndMatrices" begin include("CellScalarsVectorsAndMatricesTests.jl") end
  @time @testset "IntegrationMeshes" begin include("integrationmeshestests.jl") end
  @time @testset "Polytopes" begin include("polytopetests.jl") end
  @time @testset "RefFEs" begin include("reffetests.jl") end
  @time @testset "Meshes" begin include("meshtests.jl") end
  @time @testset "FESpaces" begin include("fespacetests.jl") end
  @time @testset "BilinearForms" begin include("bilinearformtests.jl") end
end
