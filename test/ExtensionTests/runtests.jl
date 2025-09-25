module ExtensionsTests

using Test

if !(Sys.ARCH == :aarch64 || Sys.ARCH == :i686) # Library issues on ARM and x86_64
  @testset "TikzPictures" begin include("TikzPictures.jl") end
end

end