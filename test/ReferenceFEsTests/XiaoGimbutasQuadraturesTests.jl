module XiaoGimbutasQuadraturesTests

using Test
using Gridap.ReferenceFEs, Gridap.Fields

for degree in 1:30
  quad = Quadrature(TRI,xiao_gimbutas,degree)
  @test sum(get_weights(quad)) ≈ 0.5

  ref_quad = Quadrature(TRI,duffy,degree)
  f = get_shapefuns(ReferenceFE(TRI,lagrangian,Float64,degree))
  f1 = integrate(f,get_coordinates(quad),get_weights(quad))
  f2 = integrate(f,get_coordinates(ref_quad),get_weights(ref_quad))
  err = maximum(abs.(f1 .- f2))/maximum(abs.(f1))
  # println("degree = ",degree," err = ",maximum(abs.(f1 .- f2)))
  @test err < 1.0e-7
end

for degree in 1:15
  quad = Quadrature(TET,xiao_gimbutas,degree)
  @test sum(get_weights(quad)) ≈ 0.5*1/3

  ref_quad = Quadrature(TET,duffy,degree)
  f = get_shapefuns(ReferenceFE(TET,lagrangian,Float64,degree))
  f1 = integrate(f,get_coordinates(quad),get_weights(quad))
  f2 = integrate(f,get_coordinates(ref_quad),get_weights(ref_quad))
  err = maximum(abs.(f1 .- f2))/maximum(abs.(f1))
  println("degree = ",degree," err = ",maximum(abs.(f1 .- f2)))
  @test err < 1.0e-7
end

end # module