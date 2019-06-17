module QuadraturesTests

using Test
using Gridap

const D = 2
quad = TensorProductQuadrature(orders=(2,4))

coords_ref = [
  -0.57735 0.57735 -0.57735 0.57735 -0.57735 0.57735;
  -0.774597 -0.774597 0.0 0.0 0.774597 0.774597]

@test isa(coordinates(quad),Array{Point{D,Float64},1})
@test isa(weights(quad),Array{Float64,1})

coords = coordinates(quad)

coords_arr = reshape( reinterpret(Float64,coords), (D,length(coords)) )

@test isapprox(coords_ref,coords_arr,rtol=10-3)

weigs = weights(quad)

weigs_ref = [0.555556, 0.555556, 0.888889, 0.888889, 0.555556, 0.555556]

@test isapprox(weigs,weigs_ref,rtol=10-3)

quad = Quadrature((HEX_AXIS,HEX_AXIS,HEX_AXIS),order=2)

@test isa(quad,TensorProductQuadrature{3})

end # module QuadraturesTests
