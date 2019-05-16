module Testers

using Test

using Gridap
using Gridap.Maps
using Gridap.CellMaps
using Gridap.CellValues
using Gridap.CellValues.Testers

export test_cell_map_without_gradient
export test_cell_map_with_gradient

function test_cell_map_without_gradient(
  m::CellMap{S,M,T,N},
  a::CellArray{S,M},
  b::AbstractArray{Array{T,N}} ) where {S,M,T,N}

  @test length(m) == length(b)
  @test length(m) == length(a)

  c = evaluate(m,a)
  test_iter_cell_array(c,b)

  for (mi,ai,bi) in zip(m,a,b)
    @assert isa(mi,Map{S,M,T,N})
    @assert evaluate(mi,ai) == bi
  end

  rs = maximum([size(bi) for bi in b])
  @test rs == return_size(m,cellsize(a))

end

function test_cell_map_with_gradient(
  m::CellMap{S,M,T,N},
  a::CellArray{S,M},
  b::AbstractArray{Array{T,N}},
  c::AbstractArray{Array{G,N}} ) where {S,M,T,N,G}

  test_cell_map_without_gradient(m,a,b)

  g = gradient(m)

  test_cell_map_without_gradient(g,a,c)

end

end # module Testers
