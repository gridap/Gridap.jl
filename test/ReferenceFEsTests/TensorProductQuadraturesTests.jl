module TensorProductQuadraturesTests

using Test
using Gridap.Fields, Gridap.Polynomials, Gridap.ReferenceFEs

degrees = (2, 2)
quad = Quadrature(QUAD, tensor_product, degrees)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

quad = Quadrature(QUAD, tensor_product, degrees, T=Float32)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

D = 2
degree = 3
quad = Quadrature(QUAD, tensor_product, degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

quad = Quadrature(QUAD, tensor_product, degree, T=Float32)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

@test isa(get_name(quad), String)

_exact_integral_ncube(I) = inv(prod(I))

for p in (QUAD, HEX)
  for degree in 1:10
    quad = Quadrature(p, tensor_product, degree)
    f = MonomialBasis{num_dims(p)}(Float64, degree)
    Qf = integrate(f, get_coordinates(quad), get_weights(quad))
    Ef = [_exact_integral_ncube(Tuple(term)) for term in f.terms]
    err = maximum(abs.(Qf .- Ef)) / maximum(abs.(Ef))
    @test err < 1.0e-15
  end
end

end # module
