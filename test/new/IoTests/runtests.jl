module IoTests

using Test

@testset "Interfaces" begin include("IoInterfacesTests.jl") end

@testset "Json" begin include("JsonTests.jl") end

end # module
