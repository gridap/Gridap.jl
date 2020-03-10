module RestrictedDiscreteModelsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Test

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
cell_to_mask = collect1d(apply(is_in,oldcell_to_coods))
cell_to_oldcell = findall(cell_to_mask)

model = RestrictedDiscreteModel(oldmodel,cell_to_oldcell)
test_discrete_model(model)

trian = get_triangulation(model)
@test isa(trian,RestrictedTriangulation)

trian = Triangulation(model)
@test isa(trian,RestrictedTriangulation)

model = DiscreteModel(oldmodel,cell_to_oldcell)
test_discrete_model(model)
@test isa(model,RestrictedDiscreteModel)

model = DiscreteModel(oldmodel,cell_to_mask)
test_discrete_model(model)
@test isa(model,RestrictedDiscreteModel)

#using Gridap.Visualization
#
#writevtk(model,"model")
#writevtk(trian,"trian")

end # module
