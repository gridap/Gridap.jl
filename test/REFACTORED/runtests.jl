module GridapTests

using Test
using Gridap

@time @testset "Utils" begin include("UtilsTests/runtests.jl") end

@time @testset "CellValues" begin include("CellValuesTests/runtests.jl") end

@time @testset "Fields" begin include("FieldsTests/runtests.jl") end

@time @testset "RefFEs" begin include("RefFEsTests/runtests.jl") end

@time @testset "Integration" begin include("IntegrationTests/runtests.jl") end

end # module
