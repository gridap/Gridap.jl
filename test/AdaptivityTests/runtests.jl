module RefinementTests

using Test

@testset "AdaptedGeometry" begin
  include("AdaptedGeometryTests.jl")
end

@testset "GridTransfer" begin
  include("GridTransferTests.jl")
end

end # module