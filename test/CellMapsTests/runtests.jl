module CellMapsTests

using Test

@testset "Operations" begin
  include("OperationsTests.jl")
end

@testset "ConstantCellMaps" begin
  include("ConstantCellMapsTests.jl")
end

end # module CellMapsTests
