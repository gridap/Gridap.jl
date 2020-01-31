module GridPortionsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry


domain = (-1,1,-1,1)
partition = (10,10)
oldgrid = CartesianGrid(domain,partition)

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

oldcell_to_coods = get_cell_coordinates(oldgrid)

cell_to_oldcell = findall(collect1d(apply(is_in,oldcell_to_coods)))

grid = GridPortion(oldgrid,cell_to_oldcell)
test_grid(grid)

#using Gridap.Visualization
#
#write_vtk_file(
#  grid,
#  "grid",
#  celldata=["oldcell"=>grid.cell_to_oldcell],
#  nodaldata=["oldnode"=>grid.node_to_oldnode])

end # module
