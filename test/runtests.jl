using Numa
using Test

@testset "Numa.jl" begin
  @time @testset "FieldValues" begin include("FieldValuesTests.jl") end
  @time @testset "Quadratures" begin include("QuadraturesTests.jl") end
  @time @testset "Polynomials" begin include("PolynomialsTests.jl") end
  @time @testset "CellArrays" begin include("CellArraysTests.jl") end
  @time @testset "CellQuadratures" begin include("CellQuadraturesTests.jl") end
  @time @testset "CellFunctions" begin include("CellFunctionsTests.jl") end
  #@time @testset "CellBases" begin include("CellBasesTests.jl") end
  #@time @testset "CellFields" begin include("CellFieldsTests.jl") end
  #@time @testset "CellScalarsVectorsAndMatrices" begin include("CellScalarsVectorsAndMatricesTests.jl") end
  #@time @testset "EvaluableCellArrays" begin include("EvaluableCellArraysTests.jl") end
  @time @testset "IntegrationMeshes" begin include("IntegrationMeshesTests.jl") end
  @time @testset "Polytopes" begin include("PolytopesTests.jl") end
  @time @testset "RefFEs" begin include("RefFEsTests.jl") end
  @time @testset "Meshes" begin include("MeshesTests.jl") end
  @time @testset "FESpaces" begin include("FESpacesTests.jl") end
  # @time @testset "BilinearForms" begin include("BilinearFormsTests.jl") end
end
