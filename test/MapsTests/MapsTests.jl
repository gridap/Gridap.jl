module MapsTests
##
using Gridap
using Test
using Gridap.Maps
using Gridap.Maps.Testers
using Gridap.FieldValues

import Gridap.FieldValues: inner, outer

# Unary ops

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

# Binary ops

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

# Compose

f(p::Point{2}) = 2*p
gradf(p::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
map = MockMap(a)
umap = compose(f,map)
res = evaluate(map,p)
ao = [f(ri) for ri in res]
go = [gradf(pj) for pj in p]
test_map_with_gradient(umap,p,ao,go)

map = MockMap(a)
geomap = MockMap(b)
f(p::Point{2},u::Point{2}) = 2*p + 3*u
gradf(p::Point{2},u::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
cemap = compose(f,geomap,map)
x = evaluate(geomap,p)
res = evaluate(map,p)
ao = f.(x,res)
go = gradf.(p)
test_map_with_gradient(cemap,p,ao,go)

# varinner

include("MockBasis.jl")

bas = MockBasis(a,3)
ao = [  p[i]*j+a  for j in 1:3, i in 1:length(p)]
go = [  p[i]*j  for j in 1:3, i in 1:length(p)]
test_map_with_gradient(bas,p,ao,go)

fie = MockMap(a)
scal = varinner(fie,fie)
ao = [  inner(a+pj,a+pj) for pj in p  ]
test_map_without_gradient(scal,p,ao)

bas = MockBasis(a,3)
vec = varinner(bas,fie)
ao = [  inner(p[i]*j+a,a+p[i])  for j in 1:3, i in 1:length(p)]
test_map_without_gradient(vec,p,ao)

mat = varinner(bas,bas)
ao = [ inner(p[i]*k+a,p[i]*j+a) for k in 1:3, j in 1:3, i in 1:length(p) ]
test_map_without_gradient(mat,p,ao)

# lincomb

coefs = [1.0,1.0,1.0]
ffe = lincomb(bas,coefs)
r1 = evaluate(ffe.basis,p)
ao = [ sum(r1[:,i]) for i in 1:length(p)]
test_map_without_gradient(ffe,p,ao)

# AnalyticalField

f(p::Point{2}) = 2*p
gradf(p::Point{2}) = VectorValue(2.0,2.0)
gradient(::typeof(f)) = gradf
map = AnalyticalField(f,2)
ao = f.(p)
go = gradf.(p)
test_map_with_gradient(map,p,ao,go)

# attachgeomap

geomap = MockGeoMap(2)
ao = 2*p
go = [ TensorValue(2.0,0.0,0.0,2.0) for i in 1:length(p) ]
test_map_with_gradient(geomap,p,ao,go)

bas = MockBasis(a,3)
geomap = MockGeoMap(2)
map = attachgeomap(bas,geomap)
gradbas = gradient(bas)
jaco = gradient(geomap)
rs = evaluate(bas,p)
grs = evaluate(gradient(bas),p)
js = evaluate(jaco,p)
gs = Matrix{Point{2}}(undef,(3,length(p)))
ndofs, npoints = size(grs)
for j in 1:npoints
  for i in 1:ndofs
  gs[i,j] = inv(js[j])*grs[i,j]
  end
end
test_map_with_gradient(map,p,rs,gs)

# More tests with scalar, vector, and tensor values

s_s = 5.0
s_v = VectorValue(1.1,2.3)
s_t = TensorValue(1.0,2.0,0.0,1.0)
v_s = fill(s_s,3)
v_v = fill(s_v,3)
v_t = fill(s_t,3)
a_s = fill(s_s,(3,5))
a_v = fill(s_v,(3,5))
a_t = fill(s_t,(3,5))
m_vs = TestMap(v_s,2)
m_vv = TestMap(v_v,2)
m_vt = TestMap(v_t,2)
m_as = TestMap(a_s,2)
m_av = TestMap(a_v,2)
m_at = TestMap(a_t,2)

s2_s = 1.0
s2_v = VectorValue(5.1,2.6)
s2_t = TensorValue(1.0,2.1,0.0,2.0)
v2_s = fill(s2_s,3)
v2_v = fill(s2_v,3)
v2_t = fill(s2_t,3)
a2_s = fill(s2_s,(3,5))
a2_v = fill(s2_v,(3,5))
a2_t = fill(s2_t,(3,5))
m2_vs = TestMap(v2_s,2)
m2_vv = TestMap(v2_v,2)
m2_vt = TestMap(v2_t,2)
m2_as = TestMap(a2_s,2)
m2_av = TestMap(a2_v,2)
m2_at = TestMap(a2_t,2)

test_map_without_gradient(m_av,p,a_v)

for op in (:+, :-)
  @eval begin
    for mi in [m_vs, m_vv, m_vt]
      umap = $op(mi)
      test_map_without_gradient(umap,p,$op.(mi.val))
    end
  end
end

using Gridap.CellValues.Operations: _custom_broadcast

pairs = [
  (m_vs,m2_vs), (m_vs,m2_vv), (m_vs,m2_vt),
  (m2_vs,m_vs), (m2_vv,m_vs), (m2_vt,m_vs),
  (m_vt,m2_vs), (m_vt,m2_vt)]

for op in (:+,:-,:*)
  @eval begin
    for (mi,mj) in pairs
      umap = $op(mi,mj)
      r = _custom_broadcast($op,mi.val,mj.val)
      test_map_without_gradient(umap,p,r)
    end
  end
end

pairs = [(m_vs,m2_vs),(m_vv,m2_vv),(m_vt,m2_vt)]

for op in (:(inner),)
  @eval begin
    for (mi,mj) in pairs
      umap = $op(mi,mj)
      r = _custom_broadcast($op,mi.val,mj.val)
      test_map_without_gradient(umap,p,r)
    end
  end
end

function _lincomb(a::Matrix{T},b::Vector{S}) where{T,S}
  R = Base._return_type(outer,Tuple{T,S})
  ndofs, npoints = size(a)
  v = zeros(R,npoints)
  for j in 1:npoints
    for i in 1:ndofs
      v[j] += outer(a[i,j],b[i])
    end
  end
  v
end

pairs = [
  (m_as,v2_s),(m_as,v2_v),(m_as,v2_t),
  (m_av,v2_s),(m_at,v2_s)]

for (mi,mj) in pairs
  umap = lincomb(mi,mj)
  r = _lincomb(mi.val,mj)
  test_map_without_gradient(umap,p,r)
end

end # module Maps
