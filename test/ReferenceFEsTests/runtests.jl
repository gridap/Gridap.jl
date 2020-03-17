module ReferenceFEsTests

using Test

@testset "Polytopes" begin include("PolytopesTests.jl") end

@testset "ExtrusionPolytopes" begin include("ExtrusionPolytopesTests.jl") end

@testset "Dofs" begin include("DofsTests.jl") end

@testset "MockDofs" begin include("MockDofsTests.jl") end

@testset "LagrangianDofBases" begin include("LagrangianDofBasesTests.jl") end

@testset "ReferenceFEInterfaces" begin include("ReferenceFEInterfacesTests.jl") end

@testset "NodalReferenceFEs" begin include("NodalReferenceFEsTests.jl") end

@testset "LagrangianRefFEs" begin include("LagrangianRefFEsTests.jl") end

@testset "SerendipityRefFEs" begin include("SerendipityRefFEsTests.jl") end

@testset "PDiscRefFEs" begin include("PDiscRefFEsTests.jl") end

@testset "RaviartThomasRefFEs" begin include("RaviartThomasRefFEsTests.jl") end

@testset "NedelecRefFEs" begin include("NedelecRefFEsTests.jl") end

end # module
