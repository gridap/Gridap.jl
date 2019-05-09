module MapsTests
##
using Numa
using Test
using Numa.Maps
using Numa.Maps.Testers
using Numa.FieldValues

import Numa: evaluate, gradient
import Numa: evaluate!, return_size
import Base: +, -, *, /, ∘
import Numa.FieldValues: inner, outer

include("MockMap.jl")

a = Point{2}(10,10)
b = Point{2}(15,20)
p1 = Point{2}(1,1)
p2 = Point{2}(2,2)
p3 = Point{2}(3,3)
p = [p1,p2,p3]
ao = [  a+pj for pj in p  ]
map = MockMap(a)
test_map_with_gradient(map,p,ao,p)

res = evaluate(map,p)
for op in (:+, :-)
  @eval begin
    umap = $op(map)
    ao = [ $op(rj) for rj in res]
    test_map_without_gradient(umap,p,ao)
  end
end

map1 = MockMap(a)
map2 = MockMap(b)
res1 = evaluate(map1,p)
res2 = evaluate(map2,p)
for op in (:+, :-, :inner, :outer)
  @eval begin
    umap = $op(map1,map2)
    ao = [ $op(r1i,r2i) for (r1i,r2i) in zip(res1,res2) ]
    test_map_without_gradient(umap,p,ao)
  end
end

using Numa.Maps: FieldFromCompose

f(p::Point{2}) = 2*p
gradf(p::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
map = MockMap(a)
umap = FieldFromCompose(f,map)
res = evaluate(map,p)
ao = [f(ri) for ri in res]
go = [gradf(pj) for pj in p]
test_map_with_gradient(umap,p,ao,go)

using Numa.Maps: Geomap
using Numa.Maps: FieldFromComposeExtended

map = MockMap(a)
geomap = MockMap(b)
f(p::Point{2},u::Point{2}) = 2*p + 3*u
gradf(p::Point{2},u::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
cemap = FieldFromComposeExtended(f,geomap,map)
x = evaluate(geomap,p)
res = evaluate(map,x)
ao = f.(x,res)
go = gradf.(p)
test_map_with_gradient(cemap,p,ao,go)

include("MockBasis.jl")

bas = MockBasis(a,3)
ao = [  p[i]*j+a  for j in 1:3, i in 1:length(p)]
go = [  p[i]*j  for j in 1:3, i in 1:length(p)]
test_map_with_gradient(bas,p,ao,go)

using Numa.Maps: FieldFromExpand

coefs = [1.0,1.0,1.0]
ffe = FieldFromExpand(bas,coefs)
r1 = evaluate(ffe.basis,p)
ao = [ sum(r1[:,i]) for i in 1:length(p)]
test_map_without_gradient(ffe,p,ao)

f(p::Point{2}) = 2*p
gradf(p::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
map = AnalyticalField(f,2)
ao = f.(p)
go = gradf.(p)
test_map_with_gradient(map,p,ao,go)

end # module Maps
