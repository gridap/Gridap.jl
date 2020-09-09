module CellQuadraturesTests

using Test
using FillArrays
using Gridap.Helpers
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Integration
using Gridap.CellData

using Gridap.Geometry: GridMock

degree = 3
quad = Quadrature(HEX,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

degree = (1,2,3)
quad = Quadrature(HEX,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

degree = 3
quad = Quadrature(TET,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 0.5*1/3

trian = GridMock()

degree = 3
quad = reference_cell_quadrature(trian,degree)
q = get_coordinates(quad)
@test isa(get_array(q),CompressedArray)
w = get_weights(quad)
@test isa(w,CompressedArray)

q2x = get_cell_map(trian)
x = evaluate(q2x,q)

#using Gridap.Visualization
#writevtk(trian,"trian")
#writevtk(x,"x",nodaldata=["w"=>w])

domain = (0,1,0,1)
partition = (2,3)
trian = CartesianGrid(domain,partition)
quad = reference_cell_quadrature(trian,degree)
q = get_coordinates(quad)
@test isa(get_array(q),Fill) || isa(get_array(q),CompressedArray)
w = get_weights(quad)
@test isa(w,Fill) || isa(w,CompressedArray)

dΩ = CellQuadrature(trian,degree)

vol = sum(∫(1)*dΩ)
@test vol ≈ 1

domain = (0,2,0,2,0,2)
partition = (2,3,2)
trian = CartesianGrid(domain,partition)
dΩ = CellQuadrature(trian,degree)

vol = sum(∫(1)*dΩ)
@test vol ≈ 2^3

@test isa(get_array(quad),AbstractArray{<:Quadrature})

end # module
