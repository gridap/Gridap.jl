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
@test contains(quad_tensorprod.name, "Tensor product")
@test eltype(quad_tensorprod.coordinates) == Point{num_dims(QUAD),Float32}
@test eltype(quad_tensorprod.weights) == Float32

tri_wv = Quadrature(TRI,2,T=Float32)
@test contains(tri_wv.name, "Witherden-Vincent")
@test eltype(tri_wv.coordinates) == Point{num_dims(TRI),Float32}
@test eltype(tri_wv.weights) == Float32

tet_wv = Quadrature(TET,2,T=Float32)
@test contains(tet_wv.name, "Witherden-Vincent")
@test eltype(tet_wv.coordinates) == Point{num_dims(TET),Float32}
@test eltype(tet_wv.weights) == Float32

tet_xg = Quadrature(TET,11,T=Float32)
@test contains(tet_xg.name, "Xiao-Gimbutas")

tet_duffy = Quadrature(TET,16,T=Float32)
@test contains(tet_duffy.name, "Duffy")

SIMPL4 = ExtrusionPolytope(tfill(TET_AXIS,Val(4)))
simpl4_duffy = Quadrature(SIMPL4,2,T=Float32)
@test contains(simpl4_duffy.name, "Duffy")
@test eltype(simpl4_duffy.coordinates) == Point{num_dims(SIMPL4),Float32}
@test eltype(simpl4_duffy.weights) == Float32

end # module
