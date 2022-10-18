module QuadraturesTests

using Gridap.Fields
using Gridap.ReferenceFEs
using FillArrays

coords = Point{2,Float64}[(0.0,0.0),(0.5,0.5),(1.0,1.0)]
weights = [0.5, 1.0, 0.5]
quad = GenericQuadrature(coords,weights)
test_quadrature(quad)

coords = Point{2,Float64}[(0.0,0.0),(0.5,0.5),(1.0,1.0)]
weights = Fill(1.0,3)
quad = GenericQuadrature(coords,weights)
test_quadrature(quad)

end # module
