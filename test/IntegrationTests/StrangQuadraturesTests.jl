module StrangQuadraturesTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Integration

for k in keys(Integration._strang_tri_k2n)
  reffe = ReferenceFE(TRI,lagrangian,Float64,k)
  quad1 = Quadrature(TRI,strang,k)
  quad2 = Quadrature(TRI,duffy,k)
  f = get_shapefuns(reffe)
  f1 = integrate(f,get_coordinates(quad1),get_weights(quad1))
  f2 = integrate(f,get_coordinates(quad2),get_weights(quad2))
  @test maximum(f1 .- f2)/maximum(f1) < 1.0e-7
end

for k in keys(Integration._strang_tet_k2n)
  reffe = ReferenceFE(TET,lagrangian,Float64,k)
  quad1 = Quadrature(TET,strang,k)
  quad2 = Quadrature(TET,duffy,k)
  f = get_shapefuns(reffe)
  f1 = integrate(f,get_coordinates(quad1),get_weights(quad1))
  f2 = integrate(f,get_coordinates(quad2),get_weights(quad2))
  @test maximum(f1 .- f2)/maximum(f1) < 1.0e-7
end

end # module
