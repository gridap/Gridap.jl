module ODEToolsTests

using Test

@testset "DiffOperators" begin include("DiffOperatorsTests.jl") end

@testset "ODEOperators" begin include("ODEOperatorsTests.jl") end

@testset "ODESolvers" begin include("ODESolversTests.jl") end

end # module
