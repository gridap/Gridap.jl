using Numa
using Test

@testset "Numa.jl" begin
  @time @testset "Polynomial" begin include("polynomialtests.jl") end
  #@time @testset "Polytope" begin include("polytopetests.jl") end
  #@time @testset "RefFE" begin include("reffetests.jl") end
  #@time @testset "Mesh" begin include("meshtests.jl") end
end
