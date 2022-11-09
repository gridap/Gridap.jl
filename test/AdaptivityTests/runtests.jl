module RefinementTests

using Test

@testset "AdaptedGeometry" begin
  include("AdaptedGeometryTests.jl")
end

@testset "ChangeDomain" begin
  include("ChangeDomainTests.jl")
end

end # module