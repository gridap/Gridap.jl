module Testers

using Test
using Gridap
using Gridap.Maps

export test_map_without_gradient
export test_map_with_gradient

function test_map_without_gradient(
  m::Map{S,M,T,N}, a::Array{S,M}, b::Array{T,N}) where {S,M,T,N}

  @test return_size(m,size(a)) == size(b)

  c = evaluate(m,a)

  @test c == b

  evaluate!(m,a,c)

  @test c == b

end

function test_map_with_gradient(
  m::Map{S,M,T,N}, a::Array{S,M}, b::Array{T,N}, c::Array{G,N}) where {S,M,T,N,G}

  test_map_without_gradient(m,a,b)

  g = gradient(m)

  test_map_without_gradient(g,a,c)

end

end # module Testers
