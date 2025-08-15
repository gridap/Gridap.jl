module VisualizationTests

using Test

@testset "Vtk" begin include("VtkTests.jl") end
@testset "TikzPictures" begin include("TikzPictures.jl") end

end # module
