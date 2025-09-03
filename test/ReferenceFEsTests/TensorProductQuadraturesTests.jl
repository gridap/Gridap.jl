module TensorProductQuadraturesTests

using Gridap.ReferenceFEs, Gridap.TensorValues

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

function rational_quad(degree)
  _x, w = ReferenceFEs.rational_gauss_legendre_quadrature(degree)
  x = map(VectorValue{1,eltype(_x)},_x)
  GenericQuadrature(x,w,"Rational Gauss-Legendre quadrature of degree $degree")
end
quads_1d = [rational_quad(degree) for degree in [4,3]]
quad = Quadrature(QUAD,tensor_product,quads_1d)

end # module
