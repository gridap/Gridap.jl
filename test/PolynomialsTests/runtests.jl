module PolynomialsTests

using Test

@testset "MonomialBases" begin include("MonomialBasesTests.jl") end

@testset "QGradMonomialBases" begin include("QGradMonomialBasesTests.jl") end

@testset "QCurlGradMonomialBases" begin include("QCurlGradMonomialBasesTests.jl") end

@testset "PCurlGradMonomialBases" begin include("PCurlGradMonomialBasesTests.jl") end

@testset "ModalC0Bases" begin include("ModalC0BasesTests.jl") end

@testset "JacobiPolynomialBases" begin include("JacobiPolynomialBasesTests.jl") end

#@testset "ChangeBasis" begin include("ChangeBasisTests.jl") end

end # module
