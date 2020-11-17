module AppendedTriangulationsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.Integration
using LinearAlgebra: ⋅
#using Gridap.CellData

domain = (0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)

ncells = num_cells(model)
nin = ceil(Int,2*ncells/3)
cell_to_mask = fill(false,ncells)
cell_to_mask[1:nin] .= true

grid = get_grid(model)

trian_in = RestrictedTriangulation(grid,cell_to_mask)
test_triangulation(trian_in)

trian_out = RestrictedTriangulation(grid,.! cell_to_mask)
test_triangulation(trian_out)

trian = lazy_append(trian_out,trian_in)
test_triangulation(trian)

@test isa(get_cell_coordinates(trian),AppendedArray)

@test isa(get_cell_ref_coordinates(trian),AppendedArray)

@test isa(get_cell_map(trian),AppendedArray)

@test isa(get_cell_ref_map(trian),AppendedArray)

@test !isa(get_cell_reffe(trian),AppendedArray)

@test isa(get_cell_shapefuns(trian),AppendedArray)

@test isa(get_cell_type(trian),AppendedArray)

@test isa(get_cell_id(trian),AppendedArray)

#order = 1
#quad = CellQuadrature(trian,2*order)
#quad_in = CellQuadrature(trian_in,2*order)
#quad_out = CellQuadrature(trian_out,2*order)
#
#q = get_coordinates(quad)
#w = get_weights(quad)
#@test isa(q,AppendedArray)
#@test isa(w,AppendedArray)
#
## Append triangulations of different cell type
#
#domain = (0,1,0,1)
#partition = (10,10)
#grid1 = CartesianGrid(domain,partition)
#
#domain = (1,2,0,1)
#partition = (10,10)
#grid2 = simplexify(CartesianGrid(domain,partition))
#
#trian = lazy_append(grid1,grid2)
#test_triangulation(trian)

end # module
