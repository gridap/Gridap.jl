module PCurlGradBasesTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.Arrays

xi = Point(4,2)
np = 1
x = fill(xi,np)

PT = Monomial

order = 2
D = 2
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradBasis(PT, Val(D),T,order)

v = V[
  (1.0,  0.0), (4.0,  0.0), (16.0, 0.0), (2.0,  0.0), (8.0,  0.0), (4.0,  0.0), # pterm ex
  (0.0,  1.0), (0.0,  4.0), (0.0, 16.0), (0.0,  2.0), (0.0,  8.0), (0.0,  4.0), # pterm ey
  (64.0,32.0), (32.0,16.0), (16.0, 8.0)] # sterms

g = G[
  (0.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0), (8.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (2.0, 4.0, 0.0, 0.0), (0.0, 4.0, 0.0, 0.0), # pterm ex
  (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (0.0, 0.0, 8.0, 0.0), (0.0, 0.0, 0.0, 1.0), (0.0, 0.0, 2.0, 4.0), (0.0, 0.0, 0.0, 4.0), # pterm ey
  (48.0, 0.0, 16.0, 16.0), (16.0, 16.0, 4.0, 16.0), (4.0, 16.0, 0.0, 12.0)] # sterms

vb = evaluate(b,x)

for (vi,vbi) in zip(v,vb)
  @test vi == vbi
end

vb = evaluate(b,xi)
@test vb == v

∇b = Broadcasting(gradient)(b)
gvb = evaluate(∇b,x)
for (vi,vbi) in zip(g,gvb)
  @test vi == vbi
end

gvb = evaluate(∇b,xi)
@test gvb == g

@test length(b) == 15
@test get_order(b) == 3

xi = Point(2,3,5)
np = 5
x = fill(xi,np)

order = 1
D = 3
T = Float64
V = VectorValue{D,T}
G = gradient_type(V,xi)
b = PCurlGradBasis(PT, Val(D),T,order)

@test length(b) == 15
@test get_order(b) == 2


# 1D

order = 0
D = 1
T = Float64
V = VectorValue{D,T}
b = PCurlGradBasis(PT,Val(D),T,order)

@test b isa CartProdPolyBasis{D,V,PT}

@test_throws AssertionError PCurlGradBasis(PT,Val(D),V,order)


end # module
