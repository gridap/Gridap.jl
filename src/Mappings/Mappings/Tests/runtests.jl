module MappingsTests

using Test

@testset "Mappings" begin include("MappingsTests.jl") end

@testset "MappedArrays" begin include("MappedArraysTests.jl") end

@testset "MappingArrays" begin include("MappingArraysTests.jl") end

end
