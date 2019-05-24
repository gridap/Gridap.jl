module RunTests

using Gridap
using Test

@time @testset "Helpers" begin include("HelpersTests.jl") end
@time @testset "FieldValues" begin include("FieldValuesTests.jl") end
@time @testset "CachedArrays" begin include("CachedArraysTests.jl") end
@time @testset "Maps" begin include("MapsTests/MapsTests.jl") end
@time @testset "Quadratures" begin include("QuadraturesTests.jl") end
@time @testset "Polynomials" begin include("PolynomialsTests.jl") end
@time @testset "CellValues" begin include("CellValuesTests/CellValuesTests.jl") end
@time @testset "CellMaps" begin include("CellMapsTests/runtests.jl") end
@time @testset "CellQuadratures" begin include("CellQuadraturesTests.jl") end
@time @testset "CellIntegration" begin include("CellIntegrationTests.jl") end
@time @testset "Polytopes" begin include("PolytopesTests.jl") end
@time @testset "Geometry" begin include("GeometryTests.jl") end
@time @testset "RefFEs" begin include("RefFEsTests.jl") end
@time @testset "Vtkio" begin include("VtkioTests.jl") end
@time @testset "Algebra" begin include("AlgebraTests/runtests.jl") end
@time @testset "FESpaces" begin include("FESpacesTests/runtests.jl") end

end # module RunTests
