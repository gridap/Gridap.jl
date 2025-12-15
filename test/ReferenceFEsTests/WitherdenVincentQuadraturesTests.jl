module WitherdenVincentQuadraturesTest

using Test
using Gridap.Fields, Gridap.Polynomials, Gridap.ReferenceFEs

testcases = (
  (TRI, duffy),
  (TET, duffy),
  (QUAD, tensor_product),
  (HEX, tensor_product)
)

# Witherden-Vincent quadratures are exact on Pk
# So check integration against monomials in Pk
filter = (e, o) -> (sum(e) <= o)

for (p, refquad) in testcases
  for degree in 1:maxdegree(p, witherden_vincent)
    quad = Quadrature(p, witherden_vincent, degree)
    ref_quad = Quadrature(p, refquad, degree)
    f = MonomialBasis{num_dims(p)}(Float64, degree, filter)
    Qf = integrate(f, get_coordinates(quad), get_weights(quad))
    Ef = integrate(f, get_coordinates(ref_quad), get_weights(ref_quad))
    err = maximum(abs.(Qf .- Ef)) / maximum(abs.(Ef))
    # println((degree, err))
    @test err < 1.0e-14
  end
end

_exact_integral_wedge((Ix, Iy, Iz)) = factorial(big(Ix - 1)) * factorial(big(Iy - 1)) / factorial(big(Ix + Iy)) / Iz
_exact_integral_pyramid((Ix, Iy, Iz)) = factorial(big(Ix + Iy)) * factorial(big(Iz)) / (factorial(big(Ix + Iy + Iz)) * Ix * Iy * Iz)

for (p, _exact_integral) in ((WEDGE, _exact_integral_wedge), (PYRAMID, _exact_integral_pyramid))
  for degree in 1:maxdegree(p, witherden_vincent)
    quad = Quadrature(p, witherden_vincent, degree)
    f = MonomialBasis{num_dims(p)}(Float64, degree, filter)
    Qf = integrate(f, get_coordinates(quad), get_weights(quad))
    Ef = [_exact_integral(Tuple(term)) for term in f.terms]
    err = maximum(abs.(Qf .- Ef)) / maximum(abs.(Ef))
    # println((degree, err))
    @test err < 3.0e-15
  end
end

end # module
