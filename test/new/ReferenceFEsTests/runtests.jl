module ReferenceFEsTests

using Test

@testset "Dofs" begin include("DofsTests.jl") end

@testset "MockDofs" begin include("MockDofsTests.jl") end

@testset "LagrangianDofBases" begin include("LagrangianDofBasesTests.jl") end

@testset "Polytopes" begin include("PolytopesTests.jl") end

@testset "ExtrusionPolytopes" begin include("ExtrusionPolytopesTests.jl") end

@testset "ReferenceFEInterfaces" begin include("ReferenceFEInterfacesTests.jl") end

@testset "LagrangianRefFEs" begin include("LagrangianRefFEsTests.jl") end

@testset "SerendipityRefFEs" begin include("SerendipityRefFEsTests.jl") end

end # module
