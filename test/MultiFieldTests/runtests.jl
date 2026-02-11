module MultiFieldTests

using Test

@testset "MultiFieldCellFields" begin include("MultiFieldCellFieldsTests.jl") end

@testset "MultiFieldFESpaces" begin include("MultiFieldFESpacesTests.jl") end

@testset "MultiFieldFEFunctions" begin include("MultiFieldFEFunctionsTests.jl") end

@testset "MultiFieldSparseMatrixAssemblers" begin include("MultiFieldSparseMatrixAssemblersTests.jl") end

@testset "MultiFieldFEOperators" begin include("MultiFieldFEOperatorsTests.jl") end

@testset "MultiFieldFESpacesWithLinearConstraints" begin include("MultiFieldFESpacesWithLinearConstraintsTests.jl") end

if Sys.WORD_SIZE == 32
  @info("Skipping MultiFieldFEAutodiff tests on 32-bit systems due to memory constraints.")
else
  @testset "MultiFieldFEAutodiff" begin include("MultiFieldFEAutodiffTests.jl") end
end

@testset "BlockSparseMatrixAssemblers" begin include("BlockSparseMatrixAssemblersTests.jl") end

@testset "StaticCondensation" begin include("StaticCondensationTests.jl") end

end # module
