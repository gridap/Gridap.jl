module RestrictedTriangulationsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry

domain = (0,1,0,1)
partition = (10,10)
oldgrid = CartesianGrid(domain,partition)

cell_to_oldcell = collect(1:34)
trian = RestrictedTriangulation(oldgrid,cell_to_oldcell)
test_triangulation(trian)

n_cells = num_cells(trian)

@test restrict(collect(1:n_cells),trian) == cell_to_oldcell
@test reindex(collect(1:n_cells),trian) == cell_to_oldcell
@test get_cell_id(trian) == cell_to_oldcell

oldcell_to_mask = fill(false,num_cells(oldgrid))
oldcell_to_mask[cell_to_oldcell] .= true

trian = RestrictedTriangulation(oldgrid,oldcell_to_mask)
test_triangulation(trian)

@test restrict(collect(1:n_cells),trian) == cell_to_oldcell
@test reindex(collect(1:n_cells),trian) == cell_to_oldcell
@test get_cell_id(trian) == cell_to_oldcell

end # module


