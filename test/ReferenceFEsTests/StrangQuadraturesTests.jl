module StrangQuadraturesTests

using Test
using Gridap.Fields, Gridap.ReferenceFEs

testcases = (
  (TRI, 1.0e-9),
  (TET, 1.0e-13)
)

for (p, tol) in testcases
  for degree in 1:maxdegree(p, strang)
    reffe = ReferenceFE(p, lagrangian, Float64, degree)
    quad = Quadrature(p, strang, degree)
    ref_quad = Quadrature(p, duffy, degree)
    f = get_shapefuns(reffe)
    Qf = integrate(f, get_coordinates(quad), get_weights(quad))
    Ef = integrate(f, get_coordinates(ref_quad), get_weights(ref_quad))
    err = maximum(abs.(Qf .- Ef)) / maximum(abs.(Ef))
    # println((degree, err))
    @test err < tol
  end
end

end # module
