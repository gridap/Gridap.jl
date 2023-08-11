module AdaptivityTests

using Test

@testset "AdaptedGeometry" begin
  include("RefinementRulesTests.jl")
  include("AdaptedGeometryTests.jl")
  include("FaceLabelingTests.jl")
end

@testset "Refinement" begin
  include("CartesianRefinementTests.jl")
  include("ComplexChangeDomainTests.jl")
  include("EdgeBasedRefinementTests.jl")
  include("FineToCoarseFieldsTests.jl")
  include("RefinementRuleBoundaryTests.jl")
end

@testset "CompositeQuadratures" begin
  include("CompositeQuadratureTests.jl")
end

@testset "MultiFields" begin
  include("MultifieldRefinementTests.jl")
end

end # module