module CellValuesTests

using Test

@testset "ConstantCellValues" begin
  include("ConstantCellValuesTests.jl")
end

@testset "WrappersTests" begin
  include("WrappersTests.jl")
end

@testset "OperationsTests" begin
  include("OperationsTests.jl")
end

#include("OperationsTests.jl")
#include("WrappersTests.jl")

end # module CellValuesTests
