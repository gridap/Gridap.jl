module GridapSubdivisionSurfacesRunTests

using Test

@time @testset "LoopSurfaces"  begin include("LoopSurfaceTests.jl") end

end # module
