module ODEsTests

using Test

@time @testset "ODETools" begin include("ODEsTests/runtests.jl") end

@time @testset "TransientFETools" begin include("TransientFEsTests/runtests.jl") end

# @time @testset "DiffEqsWrappers" begin include("DiffEqsWrappersTests/runtests.jl") end

# include("../bench/runbenchs.jl")

end #module
