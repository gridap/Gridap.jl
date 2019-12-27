module IoTests

using Tests

@time @testset "Interfaces" begin include("IoInterfacesTests.jl") end

@time @testset "Json" begin include("JsonTests.jl") end

end # module
