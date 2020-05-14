module AlgebraTests

using Test

@testset "AlgebraInterfaces" begin include("AlgebraInterfacesTests.jl") end

@testset "SparseMatrixCSC" begin include("SparseMatrixCSC.jl") end

@testset "SparseMatrixCSR" begin include("SparseMatrixCSR.jl") end

@testset "SymSparseMatrixCSR" begin include("SymSparseMatrixCSR.jl") end

@testset "NonlinearOperators" begin include("NonlinearOperatorsTests.jl") end

@testset "NonlinearSolvers" begin include("NonlinearSolversTests.jl") end

@testset "NLSolvers" begin include("NLSolversTests.jl") end

@testset "LinearSolvers" begin include("LinearSolversTests.jl") end

end # module
