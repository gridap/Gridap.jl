module DuffyQuadraturesTests

using Test
using Gridap.ReferenceFEs

degree = 1
quad = Quadrature(TRI,duffy,degree)
@test sum(get_weights(quad)) ≈ 0.5

degree = 1
quad = Quadrature(TRI,duffy,degree,T=Float32)
@test sum(get_weights(quad)) ≈ 0.5

degree = 4
quad = Quadrature(TET,duffy,degree)
@test sum(get_weights(quad)) ≈ 0.5*1/3

degree = 4
quad = Quadrature(TET,duffy,degree,T=Float32)
@test sum(get_weights(quad)) ≈ 0.5*1/3

end # module
