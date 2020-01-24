module DuffyQuadraturesTests

using Test
using Gridap.Integration

degree = 1
quad = DuffyQuadrature{2}(degree)
@test sum(get_weights(quad)) ≈ 0.5

degree = 4
quad = DuffyQuadrature{3}(degree)
@test sum(get_weights(quad)) ≈ 0.5*1/3

end # module
