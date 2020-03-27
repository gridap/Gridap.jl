module AppendedTriangulationsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Visualization
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration

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
@test isa(q,AppendedArray)
@test isa(w,AppendedArray)

V = TestFESpace(model=model,valuetype=Float64,order=order,reffe=:Lagrangian,conformity=:H1)

u(x) = x[1]+x[2]

_v = interpolate(V,u)
v = restrict(_v,trian)

e = u - v
el2 = sqrt(sum(integrate(e*e,trian,quad)))
@test el2 < 1.0e-8

_dv = get_cell_basis(V)
dv = restrict(_dv,trian)

cellmat =  integrate(∇(dv)*∇(dv),trian,quad)
@test isa(cellmat,AppendedArray)
@test isa(cellmat.a,CompressedArray)
@test isa(cellmat.b,CompressedArray)

#writevtk(trian_in,"trian_in")
#writevtk(trian_out,"trian_out")
#writevtk(trian,"trian",cellfields=["v"=>v])

# Append triangulations of different cell type

domain = (0,1,0,1)
partition = (10,10)
grid1 = CartesianGrid(domain,partition)

domain = (1,2,0,1)
partition = (10,10)
grid2 = simplexify(CartesianGrid(domain,partition))

trian = lazy_append(grid1,grid2)
test_triangulation(trian)

d = mktempdir()
f = joinpath(d,"trian")
writevtk(trian,f)
rm(d,recursive=true)

end # module
