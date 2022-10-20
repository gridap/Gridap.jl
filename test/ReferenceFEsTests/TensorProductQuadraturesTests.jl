module TensorProductQuadraturesTests

using Gridap.ReferenceFEs

using Test

degrees = (2,2)
quad = Quadrature(QUAD,tensor_product,degrees)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

quad = Quadrature(QUAD,tensor_product,degrees,T=Float32)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

D = 2
degree = 3
quad = Quadrature(QUAD,tensor_product,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

quad = Quadrature(QUAD,tensor_product,degree,T=Float32)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

@test isa(get_name(quad),String)

end # module
