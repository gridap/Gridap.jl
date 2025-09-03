module PolynomialsTests

using Test

@testset "PolynomialInterfaces" begin include("PolynomialInterfacesTests.jl") end

@testset "MonomialBases" begin include("MonomialBasesTests.jl") end

@testset "CurlConformBases" begin include("CurlConformBasesTests.jl") end

@testset "DivConformBases" begin include("DivConformBasesTests.jl") end

@testset "ModalC0Bases" begin include("ModalC0BasesTests.jl") end

@testset "LegendreBases" begin include("LegendreBasesTests.jl") end

@testset "ChebyshevBases" begin include("ChebyshevBasesTests.jl") end

@testset "BernsteinBases" begin include("BernsteinBasesTests.jl") end

@testset "BarycentricPΛBases" begin include("BarycentricPΛBases.jl") end

@testset "FEECBases" begin include("ExteriorCalculusBasesTests.jl") end

@testset "ForwardDiffTests.jl" begin include("ForwardDiffTests.jl") end

#@testset "ChangeBasis" begin include("ChangeBasisTests.jl") end

end # module
