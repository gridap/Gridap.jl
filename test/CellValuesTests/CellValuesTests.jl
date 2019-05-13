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

@testset "AppendTests" begin
  include("AppendTests.jl")
end

end # module CellValuesTests
