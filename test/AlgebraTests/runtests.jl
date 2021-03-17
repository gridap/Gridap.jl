module AlgebraTests

using Test

@testset "AlgebraInterfaces" begin include("AlgebraInterfacesTests.jl") end

@testset "SparseMatrixCSC" begin include("SparseMatrixCSCTests.jl") end

@testset "SparseMatrixCSR" begin include("SparseMatrixCSRTests.jl") end

@testset "SymSparseMatrixCSR" begin include("SymSparseMatrixCSRTests.jl") end

@testset "NonlinearOperators" begin include("NonlinearOperatorsTests.jl") end

@testset "NonlinearSolvers" begin include("NonlinearSolversTests.jl") end

@testset "NLSolvers" begin include("NLSolversTests.jl") end

@testset "LinearSolvers" begin include("LinearSolversTests.jl") end

end # module
