module QuadraturesTests

using Gridap.Fields
using Gridap.ReferenceFEs

coords = Point{2,Float64}[(0.0,0.0),(0.5,0.5),(1.0,1.0)]
weights = [0.5, 1.0, 0.5]
quad = GenericQuadrature(coords,weights)
test_quadrature(quad)

end # module
