module VisualizationTests

using Test

@testset "visualization_data" begin include("VisualizationDataTests.jl") end

@testset "PrintOpTrees" begin include("PrintOpTreesTests.jl") end

@testset "Vtk" begin include("VtkTests.jl") end

end # module
