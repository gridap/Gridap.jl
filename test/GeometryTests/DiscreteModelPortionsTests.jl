module DiscreteModelPortionsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry

domain = (0,2,0,2,0,2)
partition = (10,10,10)
oldmodel = CartesianDiscreteModel(domain,partition)

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

oldgrid = get_grid(oldmodel)
oldcell_to_coods = get_cell_coordinates(oldgrid)
oldcell_to_mask = lazy_map(is_in,oldcell_to_coods)

model = DiscreteModelPortion(oldmodel,oldcell_to_mask)
test_discrete_model(model)

#using Gridap.Visualization
#
#writevtk(model,"model")


end # module
