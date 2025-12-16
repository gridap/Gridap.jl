module DuffyQuadraturesTests

using Test
using Gridap.Fields, Gridap.Polynomials, Gridap.ReferenceFEs

degree = 1
quad = Quadrature(TRI, duffy, degree)
@test sum(get_weights(quad)) ≈ 0.5

degree = 1
quad = Quadrature(TRI, duffy, degree, T=Float32)
@test sum(get_weights(quad)) ≈ 0.5

degree = 4
quad = Quadrature(TET, duffy, degree)
@test sum(get_weights(quad)) ≈ 0.5 * 1 / 3

degree = 4
quad = Quadrature(TET, duffy, degree, T=Float32)
@test sum(get_weights(quad)) ≈ 0.5 * 1 / 3

_exact_integral_simplex(I) = prod(i -> factorial(big(i - 1)), I) / factorial(big(sum(I)))

filter = (e, o) -> (sum(e) <= o)

for p in (TRI, TET)
  for degree in 1:10
    quad = Quadrature(p, duffy, degree)
    f = MonomialBasis{num_dims(p)}(Float64, degree, filter)
    Qf = integrate(f, get_coordinates(quad), get_weights(quad))
    Ef = [_exact_integral_simplex(Tuple(term)) for term in f.terms]
    err = maximum(abs.(Qf .- Ef)) / maximum(abs.(Ef))
    # println((degree, err))
    @test err < 1.0e-15
  end
end

end # module
