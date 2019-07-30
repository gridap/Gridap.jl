module DuffyQuadraturesTests

using Test
using Gridap

D = 3
order=2
quad = DuffyQuadrature{D}(order)
@test isa(coordinates(quad),Vector{Point{D,Float64}})
@test isa(weights(quad),Vector{Float64})
@test sum(weights(quad)) ≈ 1/6

D = 2
order=2
quad = DuffyQuadrature{D}(order)
x = coordinates(quad)
r = Point{D,Float64}[
  (0.155051, 0.178559), (0.644949, 0.0750311),
  (0.155051, 0.66639), (0.644949, 0.28002)]
@test isapprox(reinterpret(Float64,x),reinterpret(Float64,r),rtol=1e-3)
@test sum(weights(quad)) ≈ 1/2

D = 1
order=2
quad = DuffyQuadrature{D}(order)
@test sum(weights(quad)) ≈ 1

end # module
