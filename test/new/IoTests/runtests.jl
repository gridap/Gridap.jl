module IoTests

using Test

@testset "Interfaces" begin include("IoInterfacesTests.jl") end

@testset "Json" begin include("JsonTests.jl") end

@testset "Bson" begin include("BsonTests.jl") end

end # module
