module VisualizationTests

using Test

@testset "Vtk" begin include("VtkTests.jl") end

@info "ARCH: $(Sys.ARCH)"
if !(Sys.ARCH === :aarch64 || Sys.ARCH === :x86_64) # Library issues on ARM and x86_64
  @testset "TikzPictures" begin include("TikzPictures.jl") end
end

end # module
