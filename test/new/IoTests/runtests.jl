module IoTests

using Test

@time @testset "Interfaces" begin include("IoInterfacesTests.jl") end

@time @testset "Json" begin include("JsonTests.jl") end

end # module
