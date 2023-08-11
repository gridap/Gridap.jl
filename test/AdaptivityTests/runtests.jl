module AdaptivityTests

using Test

@testset "Refinement" begin
  include("RefinementRulesTests.jl")
  include("AdaptedGeometryTests.jl")
  include("FaceLabelingTests.jl")
  include("CartesianRefinementTests.jl")
  include("ComplexChangeDomainTests.jl")
  include("EdgeBasedRefinementTests.jl")
  include("FineToCoarseFieldsTests.jl")
  include("RefinementRuleBoundaryTests.jl")
  include("MultifieldRefinementTests.jl")
end

@testset "CompositeQuadratures" begin
  include("CompositeQuadratureTests.jl")
end

end # module