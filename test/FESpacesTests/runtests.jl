module FESpacesTests
using Test
include("runtests_1.jl")
include("runtests_2.jl")
@testset "LinearizedFESpacesTests" begin include("LinearizedFESpacesTests.jl") end
end # module
