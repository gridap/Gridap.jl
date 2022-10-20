module QuadraturesTests

using Test
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

# tests with Float32 usage with T keyword

quad_tensorprod = Quadrature(QUAD,2,T=Float32)
@test eltype(quad_tensorprod.coordinates) == Point{num_dims(QUAD),Float32}
@test eltype(quad_tensorprod.weights) == Float32

quad_strang = Quadrature(TET,2,T=Float32)
@test eltype(quad_strang.coordinates) == Point{num_dims(TET),Float32}
@test eltype(quad_strang.weights) == Float32

quad_duffy = Quadrature(TRI,2,T=Float32)
@test eltype(quad_duffy.coordinates) == Point{num_dims(TRI),Float32}
@test eltype(quad_duffy.weights) == Float32


end # module
