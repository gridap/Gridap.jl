module AppendedTriangulationsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.Integration
using LinearAlgebra: ⋅
using Gridap.CellData
using Gridap.FESpaces

domain = (0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)

ncells = num_cells(model)
nin = ceil(Int,2*ncells/3)
cell_to_mask = fill(false,ncells)
cell_to_mask[1:nin] .= true

grid = get_grid(model)

trian_in = RestrictedTriangulation(grid,cell_to_mask)
trian_out = RestrictedTriangulation(grid,collect(Bool, .! cell_to_mask))

trian = lazy_append(trian_out,trian_in)
test_triangulation(trian)

order = 1
quad = CellQuadrature(trian,2*order)
quad_in = CellQuadrature(trian_in,2*order)
quad_out = CellQuadrature(trian_out,2*order)

q = get_coordinates(quad)
w = get_weights(quad)
@test isa(get_array(q.q),AppendedArray)
@test isa(w.refvals,AppendedArray)

V = TestFESpace(model=model,valuetype=Float64,order=order,reffe=:Lagrangian,conformity=:H1)

u(x) = x[1]+x[2]

v = interpolate(u,V)

e = u - v
el2 = sqrt(sum(integrate(e*e,quad)))
@test el2 < 1.0e-8

dv = get_cell_basis(V)

cellmat =  integrate(∇(dv)⋅∇(dv),quad)
@test isa(cellmat,AppendedArray)
@test isa(cellmat.a,CompressedArray)
@test isa(cellmat.b,CompressedArray)

end # module
