module CartesianDiscreteModelsTests

using Gridap
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: CartesianMap
using Gridap.Io
using Test

domain = (0,1,0,1,0,1)
partition = (3,3,3)
parts = nothing
model = CartesianDiscreteModel(parts,domain,partition)
test_discrete_model(model)

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)
test_discrete_model(model)
@test is_oriented(get_grid(model)) == true
@test is_oriented(get_grid_topology(model)) == true

grid_topology = get_grid_topology(model)
cell_to_pindices = get_cell_permutations(grid_topology)

labels = get_face_labeling(model)
@test num_tags(labels) == 28
@test num_entities(labels) == 27

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
model = CartesianDiscreteModel(domain,partition;map=polar)
test_discrete_model(model)
@test is_oriented(get_grid(model)) == true

model2 = from_dict(DiscreteModel,to_dict(model))
test_discrete_model(model2)

model3 = CartesianDiscreteModel(desc,CartesianIndex(2,2,2),CartesianIndex(3,3,3))
test_discrete_model(model3)

domain = (0,1,0,1)
partition = (4,4)
desc = CartesianDescriptor(domain,partition;isperiodic=(true,true))
model = CartesianDiscreteModel(desc)
@test num_vertices(model) == 16
model = CartesianDiscreteModel(desc,CartesianIndex(1,1),CartesianIndex(4,4))
@test num_vertices(model) == 16

desc = CartesianDescriptor(domain,partition;isperiodic=(true,false))
model = CartesianDiscreteModel(desc)
@test num_vertices(model) == 20
model = CartesianDiscreteModel(desc,CartesianIndex(1,1),CartesianIndex(4,4))
@test num_vertices(model) == 20

desc = CartesianDescriptor(domain,partition;isperiodic=(false,true))
model = CartesianDiscreteModel(desc)
@test num_vertices(model) == 20
model = CartesianDiscreteModel(desc,CartesianIndex(1,1),CartesianIndex(4,4))
@test num_vertices(model) == 20

desc = CartesianDescriptor(domain,partition;isperiodic=(false,false))
model = CartesianDiscreteModel(desc)
@test num_vertices(model) == 25
model = CartesianDiscreteModel(desc,CartesianIndex(1,1),CartesianIndex(4,4))
@test num_vertices(model) == 25

#using Gridap.Visualization
#writevtk(model2,"model2")

end # module
