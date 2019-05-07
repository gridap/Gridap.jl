module CellValuesTests

using Test

@testset "ConstantCellValues" begin
  include("ConstantCellValuesTests.jl")
end

#include("Mocks.jl")
#include("OperationsTests.jl")
#include("WrappersTests.jl")

end # module CellValuesTests
