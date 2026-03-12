module GridapRunTests

using Test

TESTCASE = get(ENV, "TESTCASE", "all")

if TESTCASE ∈ ("all", "unit-basics")
  @time @testset "Helpers"      begin include("HelpersTests/runtests.jl") end
  @time @testset "Io"           begin include("IoTests/runtests.jl") end
  @time @testset "Algebra"      begin include("AlgebraTests/runtests.jl") end
  @time @testset "Arrays"       begin include("ArraysTests/runtests.jl") end
  @time @testset "TensorValues" begin include("TensorValuesTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-fields")
  @time @testset "Fields"      begin include("FieldsTests/runtests.jl") end
  @time @testset "Polynomials" begin include("PolynomialsTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-referencefes")
  @time @testset "ReferenceFEs" begin include("ReferenceFEsTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-geometry")
  @time @testset "Geometry" begin include("GeometryTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-celldata")
  @time @testset "CellData" begin include("CellDataTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-visualization")
  @time @testset "Visualization" begin include("VisualizationTests/runtests.jl") end
  @time @testset "Aqua"          begin include("Aqua.jl") end
end

if TESTCASE ∈ ("all", "unit-fespaces-1")
  @time @testset "FESpaces (1/2)" begin include("FESpacesTests/runtests_1.jl") end
end

if TESTCASE ∈ ("all", "unit-fespaces-2")
  @time @testset "FESpaces (2/2)" begin include("FESpacesTests/runtests_2.jl") end
end

if TESTCASE ∈ ("all", "unit-multifield")
  @time @testset "MultiField" begin include("MultiFieldTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-odes")
  @time @testset "ODEs" begin include("ODEsTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "unit-adaptivity")
  @time @testset "Adaptivity" begin include("AdaptivityTests/runtests.jl") end
end

if TESTCASE ∈ ("all", "drivers")
  @time @testset "Drivers" begin include("GridapTests/runtests.jl") end
end

end # module
