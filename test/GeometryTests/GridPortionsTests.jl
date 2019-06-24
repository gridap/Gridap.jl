module GridPortionsTests

using Test
using Gridap

oldgrid = CartesianGrid(partition=(3,3))

n = ncells(oldgrid)
nini = ceil(n/3)
newcell_to_oldcell = collect(Int,nini:n)
grid = GridPortion(oldgrid,newcell_to_oldcell)

test_grid(grid,16,7)

#writevtk(grid,"grid")

model = CartesianDiscreteModel(partition=(3,4,2))

labels = FaceLabels(model,2)
oldgrid = Grid(model,2)

mask = [ (label == 23 || label == 24 || label == 25) for label in labels ]

newcell_to_oldcell = findall(mask)
grid = GridPortion(oldgrid,newcell_to_oldcell)

test_grid(grid,60,20)

#writevtk(grid,"grid",celldata=["oldcells"=>newcell_to_oldcell])

end # module
