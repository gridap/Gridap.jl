module UtilsTests

using Test

@testset "Helpers" begin include("HelpersTests.jl") end

@testset "CachedArrays" begin include("CachedArraysTests.jl") end

end # module
