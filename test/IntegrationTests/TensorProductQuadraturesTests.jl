module TensorProductQuadraturesTests

using Gridap.Integration

using Test

degrees = (2,2)
quad = TensorProductQuadrature(degrees)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

D = 2
degree = 3
quad = TensorProductQuadrature{D}(degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

end # module
