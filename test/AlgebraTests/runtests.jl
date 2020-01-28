module AlgebraTests

using Test

@testset "SparseMatrixCSC" begin include("SparseMatrixCSC.jl") end

@testset "SparseMatrixCSR" begin include("SparseMatrixCSR.jl") end

@testset "SymSparseMatrixCSR" begin include("SymSparseMatrixCSR.jl") end

@testset "BlockArraysCOO" begin include("BlockArraysCOOTests.jl") end

@testset "NonLinearOperators" begin include("NonLinearOperatorsTests.jl") end

@testset "NonLinearSolvers" begin include("NonLinearSolversTests.jl") end

@testset "NLSolvers" begin include("NLSolversTests.jl") end

@testset "LinearSolvers" begin include("LinearSolversTests.jl") end


end # module
