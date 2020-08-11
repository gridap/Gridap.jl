module VisualizationTests

using Test

@testset "PrintOpTrees" begin include("PrintOpTreesTests.jl") end

@testset "Vtk" begin include("VtkTests.jl") end

@testset "visualization_data" begin include("VisualizationData.jl") end

end # module
