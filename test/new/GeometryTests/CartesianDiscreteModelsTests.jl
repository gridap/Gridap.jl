module CartesianDiscreteModelsTests

using Gridap
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: CartesianMap
using Test

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)
test_discrete_model(model)
@test is_oriented(get_grid(model)) == true
@test is_oriented(get_grid_topology(model)) == true

labels = get_face_labeling(model)
@test num_tags(labels) == 29
@test num_entities(labels) == 28

map = get_cell_map(get_grid(model))
@test isa(map,CartesianMap)

desc = get_cartesian_descriptor(model)

function polar(q)
  r, t, z = q
  x = r*cos(t)
  y = r*sin(t)
  Point(x,y,z)
end

domain = (1,2,0,pi,0,0.5)
partition = (10,30,4)
model = CartesianDiscreteModel(domain,partition,polar)
test_discrete_model(model)
@test is_oriented(get_grid(model)) == true

#using Gridap.Visualization
#writevtk(model,"model")

end # module
