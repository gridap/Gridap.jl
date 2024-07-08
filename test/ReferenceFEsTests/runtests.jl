module ReferenceFEsTests

using Test

@testset "Polytopes" begin include("PolytopesTests.jl") end

@testset "ExtrusionPolytopes" begin include("ExtrusionPolytopesTests.jl") end

@testset "GeneralPolytopes" begin include("GeneralPolytopesTests.jl") end

@testset "Dofs" begin include("DofsTests.jl") end

@testset "MockDofs" begin include("MockDofsTests.jl") end

@testset "LagrangianDofBases" begin include("LagrangianDofBasesTests.jl") end

@testset "ReferenceFEInterfaces" begin include("ReferenceFEInterfacesTests.jl") end

@testset "LagrangianRefFEs" begin include("LagrangianRefFEsTests.jl") end

@testset "CLagrangianRefFEs" begin include("CLagrangianRefFEsTests.jl") end

@testset "SerendipityRefFEs" begin include("SerendipityRefFEsTests.jl") end

@testset "PDiscRefFEs" begin include("PDiscRefFEsTests.jl") end

@testset "Quadratures" begin include("QuadraturesTests.jl") end

@testset "TensorProductQuadratures" begin include("TensorProductQuadraturesTests.jl") end

@testset "DuffyQuadratures" begin include("DuffyQuadraturesTests.jl") end

@testset "StrangQuadratures" begin include("StrangQuadraturesTests.jl") end

@testset "RaviartThomasRefFEs" begin include("RaviartThomasRefFEsTests.jl") end

@testset "NedelecRefFEs" begin include("NedelecRefFEsTests.jl") end

@testset "CDLagrangianRefFEs" begin include("CDLagrangianRefFEsTests.jl") end

@testset "BezierRefFEs" begin include("BezierRefFEsTests.jl") end

@testset "ModalC0RefFEs" begin include("ModalC0RefFEsTests.jl") end

end # module
