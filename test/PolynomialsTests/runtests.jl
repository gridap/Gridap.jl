module PolynomialsTests

using Test

@testset "PolynomialInterfaces" begin include("PolynomialInterfacesTests.jl") end

@testset "MonomialBases" begin include("MonomialBasesTests.jl") end

@testset "QGradBases" begin include("QGradBasesTests.jl") end

@testset "QCurlGradBases" begin include("QCurlGradBasesTests.jl") end

@testset "PGradBases" begin include("PGradBasesTests.jl") end

@testset "PCurlGradBases" begin include("PCurlGradBasesTests.jl") end

@testset "ModalC0Bases" begin include("ModalC0BasesTests.jl") end

@testset "LegendreBases" begin include("LegendreBasesTests.jl") end

@testset "BernsteinBases" begin include("BernsteinBasesTests.jl") end

#@testset "ChangeBasis" begin include("ChangeBasisTests.jl") end

end # module
