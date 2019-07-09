module FESpacesTestsAll

using Test

@testset "FESpaces" begin include("FESpacesTests.jl") end

@testset "FEFunctions" begin include("FEFunctionsTests.jl") end

@testset "FEBases" begin include("FEBasesTests.jl") end

@testset "Assemblers" begin include("AssemblersTests.jl") end

@testset "FEterms" begin include("FETermsTests.jl") end

@testset "FEOperators" begin include("FEOperatorsTests.jl") end

@testset "NonLinearFEOperators" begin include("NonLinearFEOperatorsTests.jl") end

end # module
