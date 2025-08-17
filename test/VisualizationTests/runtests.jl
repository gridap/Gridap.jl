module VisualizationTests

using Test

@testset "Vtk" begin include("VtkTests.jl") end

if !(Sys.ARCH == :aarch64 || Sys.ARCH == :i686) # Library issues on ARM and x86_64
  @testset "TikzPictures" begin include("TikzPictures.jl") end
end

end # module
