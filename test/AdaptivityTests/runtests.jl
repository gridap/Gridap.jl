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
  include("UnstructuredUniformRefinementTests.jl")
end

@testset "CompositeQuadratures" begin
  include("CompositeQuadratureTests.jl")
end

@testset "MacroFETests" begin
  include("MacroFETests.jl")
  include("MacroFEStokesTests.jl")
end

@testset "AMR" begin
  include("AdaptiveMeshRefinementTests.jl")
end

@testset "PolytopalCoarsening" begin
  include("PolytopalCoarseningTests.jl")
end

end # module