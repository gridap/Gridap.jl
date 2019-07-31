module SimplexifyTests

using Gridap

grid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0),partition=(3,4))
tgrid = simplexify(grid)
test_grid(tgrid,20,24)
#writevtk(tgrid,"tgrid")

grid = CartesianGrid(partition=(3,4,3))
tgrid = simplexify(grid)
test_grid(tgrid,80,216)
#writevtk(tgrid,"tgrid")

model = CartesianDiscreteModel(partition=(3,4,2))
tmodel = simplexify(model)

model = CartesianDiscreteModel(partition=(3,4))
tmodel = simplexify(model)
#writevtk(tmodel,"tmodel")


end # module
