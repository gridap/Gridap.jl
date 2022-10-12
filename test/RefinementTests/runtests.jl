module RefinementTests

using Test

@testset "RefinedGeometry" begin
  include("RefinedGeometryTests.jl")
end

@testset "ChangeDomain" begin
  include("ChangeDomainTests.jl")
end

end # module