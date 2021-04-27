module StrangQuadraturesTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs

for k in keys(ReferenceFEs._strang_tri_k2n)
  reffe = ReferenceFE(TRI,lagrangian,Float64,k)
  quad1 = Quadrature(TRI,strang,k)
  quad2 = Quadrature(TRI,duffy,k)
  f = get_shapefuns(reffe)
  f1 = integrate(f,get_coordinates(quad1),get_weights(quad1))
  f2 = integrate(f,get_coordinates(quad2),get_weights(quad2))
  err = maximum(abs.(f1 .- f2))/maximum(abs.(f1))
  @test err < 1.0e-7
  #@show (2,k,err)
end

for k in keys(ReferenceFEs._strang_tet_k2n)
  reffe = ReferenceFE(TET,lagrangian,Float64,k)
  quad1 = Quadrature(TET,strang,k)
  quad2 = Quadrature(TET,duffy,k)
  f = get_shapefuns(reffe)
  f1 = integrate(f,get_coordinates(quad1),get_weights(quad1))
  f2 = integrate(f,get_coordinates(quad2),get_weights(quad2))
  err = maximum(abs.(f1 .- f2))/maximum(abs.(f1))
  @test err < 1.0e-13
  #@show (3,k,err)
end

end # module
