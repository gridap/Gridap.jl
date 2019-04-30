module MapsTests
##
using Numa
using Test
using Numa.Maps
using Numa.FieldValues

import Numa: evaluate, gradient
import Numa: evaluate!, return_size
import Base: +, -, *, /, ∘
import Numa.FieldValues: inner, outer

include("MockMap.jl")

"""
Check whether all the queries of the Map interface have been defined
"""
function is_a_map(m::Map{S,M,T,N}) where {S,M,T,N}
  sa = tuple([0 for i in 1:M]...)
  sb = return_size(m,sa)
  a = zeros(S,sa)
  b = Array{S,M}(undef,sb)
  evaluate!(m,a,b)
  #gm = gradient(m)
  true
end

a = Point{2}(10,10)
b = Point{2}(15,20)
p1 = Point{2}(1,1)
p2 = Point{2}(2,2)
p3 = Point{2}(3,3)
p = [p1,p2,p3]
##
@testset "MockMap" begin
  length(p)
  map = MockMap(a)
  @test is_a_map(map)
  res = evaluate(map,p)
  for i in 1:length(p)
    @test res[i] == a+p[i]
  end
  @test return_size(map,size(p)) == size(p)
  gmap = gradient(map)
  gres = evaluate(gmap,p)
  for i in 1:length(p)
    @test gres[i] == p[i]
  end
end

# Unary Operators
mymap = MockMap(a)
using Numa.Maps: MapFromUnaryOp
res = evaluate(mymap,p)
@testset "UnaryOp" begin
  for op in (:+, :-)
    @eval begin
      umap = $op(mymap)
      @test is_a_map(umap)
      res2 = evaluate(umap,p)
      for i in 1:length(p)
        @test res2[i] == $op(res[i])
      end
    end
  end
end

# Binary Operators
map1 = MockMap(a)
map2 = MockMap(b)
res1 = evaluate(map1,p)
res2 = evaluate(map2,p)
using Numa.Maps: MapFromBinaryOp
@testset "BinaryOp" begin
  for op in (:+, :-, :inner, :outer)
    @eval begin
      umap = MapFromBinaryOp($op,map1,map2)
      @test is_a_map(umap)
      resu = evaluate(umap,p)
      for i in 1:length(p)
        @test resu[i] == $op(res1[i],res2[i])
      end
      isa(umap,MapFromBinaryOp{typeof($op),MockMap{2}})
    end
  end
end

# Compose
f(p::Point{2}) = 2*p
gradf(p::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
@testset "ComposeField" begin
  using Numa.Maps: FieldFromCompose
  @test MockMap <: Field
  map = MockMap(a)
  umap = FieldFromCompose(f,map)
  @test is_a_map(umap)
  res = evaluate(map,p)
  resu = evaluate(umap,p)
  for i in 1:length(p)
    @test resu[i] == f(res[i])
  end
  gumap = gradient(umap)
  gres = evaluate(gumap,p)
  for i in 1:length(p)
    @test gres[i] == gradf(p[i])
  end
end

# FieldFromComposeExtended
using Numa.Maps: Geomap
MockMap <: Geomap
map = MockMap(a)
geomap = MockMap(b)
@testset "ComposeExtended" begin
  using Numa.Maps: FieldFromComposeExtended
  cemap = FieldFromComposeExtended(f,geomap,map)
  @test is_a_map(cemap)
  res = evaluate(cemap,p)
  for i in 1:length(p)
    res[i] == f(evaluate(map,evaluate(geomap,[p[i]]))...)
  end
  gcemap = gradient(cemap)
  gres = gradient(cemap)
  gres = evaluate(gcemap,p)
  for i in 1:length(p)
    @test gres[i] == gradf(p[i])
  end
end

include("MockBasis.jl")

bas = MockBasis(a,3)
@test is_a_map(bas)

@testset "MockBasis" begin
  res = evaluate(bas,p)
  v_size = return_size(bas, size(p))
  v = Array{Point{2},2}(undef, v_size)
  for j = 1:3
    for (i,pi) in enumerate(p)
      @test res[j,i] == pi*j+a
    end
  end
  gbas = gradient(bas)
  gres = evaluate(gbas,p)
  for j = 1:3
    for (i,pi) in enumerate(p)
      @test gres[j,i] == j*pi
    end
  end
end

using Numa.Maps: FieldFromExpand
@testset "FieldFromExpand" begin
  coefs = [1.0,1.0,1.0]
  ffe = FieldFromExpand(bas,coefs)
  @test is_a_map(ffe)
  res = evaluate(ffe,p)
  r1 = evaluate(ffe.basis,p)
  for i in 1:ffe.basis.dim
    @test res[i] == sum(r1[:,i])
  end
end

end # module Maps
