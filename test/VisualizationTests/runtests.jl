module VisualizationTests

using Test

@testset "visualization_data" begin include("VisualizationDataTests.jl") end

@testset "Vtk" begin include("VtkTests.jl") end

end # module
