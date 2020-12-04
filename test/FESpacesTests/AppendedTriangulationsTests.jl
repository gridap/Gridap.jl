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
using Gridap.ReferenceFEs
using FillArrays

domain = (0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)

ncells = num_cells(model)
nin = ceil(Int,2*ncells/3)
cell_to_mask = fill(false,ncells)
cell_to_mask[1:nin] .= true

grid = get_grid(model)

Ω_in = RestrictedTriangulation(grid,cell_to_mask)
Ω_out = RestrictedTriangulation(grid,.! cell_to_mask)
Ω = lazy_append(Ω_out,Ω_in)
test_triangulation(Ω)

order = 1
degree = 2*order
quad_in = CellQuadrature(Ω_in,degree)
quad_out = CellQuadrature(Ω_out,degree)
quad = CellQuadrature(Ω,degree)
@test isa(quad.trian,AppendedTriangulation)
@test isa(quad.cell_quad,AppendedArray)
@test isa(quad.cell_point,AppendedArray)
@test isa(quad.cell_weight,AppendedArray)

V = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order),conformity=:H1)

v(x) = x[1]+x[2]

vh = interpolate(v,V)

e = v - vh
el2 = sqrt(sum(integrate(e*e,quad)))
@test el2 < 1.0e-8

x = get_cell_points(quad)

dv = get_cell_shapefuns(V)
du = get_cell_shapefuns_trial(V)
@test isa(dv(x),AppendedArray)
@test isa(∇(dv)(x),AppendedArray)
@test isa(du(x),AppendedArray)
@test isa(∇(du)(x),AppendedArray)

cellmat = integrate( ∇(dv)⋅∇(du), quad )

@test isa(cellmat,AppendedArray)
@test isa(cellmat.a,Fill)
@test isa(cellmat.b,Fill)

end # module
