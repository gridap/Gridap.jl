module Grida

using Test
using Gridap

@time @testset "Utils" begin include("UtilsTests/runtests.jl") end

@time @testset "CellValues" begin include("CellValuesTests/runtests.jl") end

@time @testset "Fields" begin include("FieldsTests/runtests.jl") end

end # module
