module ODEsTests

using Test

# @time @testset "TimeDerivatives" begin include("TimeDerivativesTests.jl") end

# @time @testset "ODEOperators" begin include("ODEOperatorsTests.jl") end

# @time @testset "ODESolvers" begin include("ODESolversTests.jl") end

# @time @testset "ODESolutions" begin include("ODESolutionsTests.jl") end

# @time @testset "TransientFESpaces" begin include("TransientFESpacesTests.jl") end

# @time @testset "TransientCellFields" begin include("TransientCellFieldsTests.jl") end

# @time @testset "TransientFEOperators" begin include("TransientFEOperatorsTests.jl") end

@time @testset "TransientFESolutions" begin include("TransientFESolutionsTests.jl") end

# @time @testset "ODESolversAll" begin include("ODESolversAllTests/runtests.jl") end

# @time @testset "TransientProblems" begin include("TransientProblemsTests/runtests.jl") end

# TODO Find a way to run the same tests as in
# - TransientFETests
# - TransientFEOperatorsTests
# - AffineFEOperatorsTests
# - ConstantFEOperatorsTests
# - Transient2ndOrderFEOperatorsTests

# @time @testset "DiffEqsWrappers" begin include("_DiffEqsWrappersTests.jl") end

# include("../bench/runbenchs.jl")

end # module ODEsTests
