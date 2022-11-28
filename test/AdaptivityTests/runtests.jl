module AdaptivityTests

using Test

@testset "AdaptedGeometry" begin
  include("AdaptedGeometryTests.jl")
end

@testset "GridTransfer" begin
  include("GridTransferTests.jl")
end

@testset "CompositeQuadratures" begin
  include("CompositeQuadratureTests.jl")
end

@testset "MultiFields" begin
  include("MultifieldGridTransferTests.jl")
end

end # module