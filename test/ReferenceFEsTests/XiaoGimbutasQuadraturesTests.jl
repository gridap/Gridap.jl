module XiaoGimbutasQuadraturesTests

using Test
using Gridap.Fields, Gridap.Polynomials, Gridap.ReferenceFEs

testcases = (
  (TRI, duffy),
  (TET, duffy),
  (QUAD, tensor_product),
  (HEX, tensor_product)
)

# Xiao-Gimbutas quadratures are exact on Pk
# So check integration against monomials in Pk
filter = (e, o) -> (sum(e) <= o)

for (p, refquad) in testcases
  for degree in 1:maxdegree(p, xiao_gimbutas)
    quad = Quadrature(p, xiao_gimbutas, degree)
    ref_quad = Quadrature(p, refquad, degree)
    f = MonomialBasis{num_dims(p)}(Float64, degree, filter)
    Qf = integrate(f, get_coordinates(quad), get_weights(quad))
    Ef = integrate(f, get_coordinates(ref_quad), get_weights(ref_quad))
    err = maximum(abs.(Qf .- Ef)) / maximum(abs.(Ef))
    # println((degree, err))
    @test err < 1.0e-14
  end
end

end # module
