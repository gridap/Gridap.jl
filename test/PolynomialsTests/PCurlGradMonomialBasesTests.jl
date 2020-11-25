module PCurlGradMonomialBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.Arrays

xi = Point(4,2)
np = 1
x = fill(xi,np)

order = 2
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradMonomialBasis{D}(T,order)

v = V[
  (1.0, 0.0), (0.0, 1.0), (4.0, 0.0), (0.0, 2.0), (16.0, 0.0), (0.0, 4.0),
  (2.0, 0.0), (0.0, 4.0), (8.0, 0.0), (0.0, 8.0), (4.0, 0.0), (0.0, 16.0),
  (64.0, 32.0), (32.0, 16.0), (16.0, 8.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0),
  (0.0, 0.0, 0.0, 1.0), (8.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 4.0),
  (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (2.0, 4.0, 0.0, 0.0),
  (0.0, 0.0, 2.0, 4.0), (0.0, 4.0, 0.0, 0.0), (0.0, 0.0, 8.0, 0.0),
  (48.0, 0.0, 16.0, 16.0), (16.0, 16.0, 4.0, 16.0), (4.0, 16.0, 0.0, 12.0)]

  vb = evaluate(b,x)

  for (vi,vbi) in zip(v,vb)
    @test vi == vbi
  end

  ∇b = Broadcasting(gradient)(b)
  gvb = evaluate(∇b,x)
  for (vi,vbi) in zip(g,gvb)
    @test vi == vbi
  end

  @test num_terms(b) == 15
  @test get_order(b) == 2

  xi = Point(2,3,5)
  np = 5
  x = fill(xi,np)

  order = 1
  D = 3
  T = Float64
  V = VectorValue{D,T}
  G = gradient_type(V,xi)
  b = PCurlGradMonomialBasis{D}(T,order)

  @test num_terms(b) == 15
  @test get_order(b) == 1

end # module
