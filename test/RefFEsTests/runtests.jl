module RefFEsTestsAll

using Test

@testset "PolynomialBases" begin include("PolynomialBasesTests.jl") end

@testset "Polytopes" begin include("PolytopesTests.jl") end

@testset "DOFBases" begin include("DOFBasesTests.jl") end

@testset "RefFEs" begin include("RefFEsTests.jl") end

@testset "PDiscRefFEs" begin include("PDiscRefFEsTests.jl") end

@testset "GenericRefFEsTests" begin include("GenericRefFEsTests.jl") end

end # module
