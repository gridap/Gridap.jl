module GridapSubdivisionSurfacesRunTests

using Test

@time @testset "TorusLoopSurfaces"  begin include("TorusLoopSurfaceTests.jl") end
@time @testset "SphereLoopSurfaces"  begin include("SphereLoopSurfaceTests.jl") end

end # module
