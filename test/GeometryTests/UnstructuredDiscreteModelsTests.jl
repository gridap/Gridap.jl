module UnstructuredDiscreteModelsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Io

m = DiscreteModelMock()
g = get_grid(m)

model = UnstructuredDiscreteModel(g)
test_discrete_model(model)

domain = (0,1,0,1)
partition = (2,3)
m = CartesianDiscreteModel(domain,partition)
model = UnstructuredDiscreteModel(m)
test_discrete_model(model)
@test model === UnstructuredDiscreteModel(model)

m = DiscreteModelMock()
model = UnstructuredDiscreteModel(m)
test_discrete_model(model)
@test model === UnstructuredDiscreteModel(model)

m = DiscreteModelMock()
model = UnstructuredDiscreteModel(m)

dict = to_dict(model)
model2 = from_dict(UnstructuredDiscreteModel,dict)
test_discrete_model(model2)

model2 = from_json(UnstructuredDiscreteModel,to_json(model))
test_discrete_model(model2)

d = mktempdir()
f = joinpath(d,"model.jld2")

to_jld2_file(model,f)
model2 == from_jld2_file(typeof(model),f)

@test model.face_labeling.d_to_dface_to_entity == model2.face_labeling.d_to_dface_to_entity
@test model.face_labeling.tag_to_entities == model2.face_labeling.tag_to_entities
@test model.face_labeling.tag_to_name == model2.face_labeling.tag_to_name
@test model.grid.cell_node_ids == model2.grid.cell_node_ids
@test model.grid.cell_types == model2.grid.cell_types
@test model.grid.node_coordinates == model2.grid.node_coordinates
@test model.grid.reffes == model2.grid.reffes
@test model.grid_topology.cell_type == model2.grid_topology.cell_type
@test model.grid_topology.polytopes == model2.grid_topology.polytopes
@test model.grid_topology.vertex_coordinates == model2.grid_topology.vertex_coordinates

rm(d,recursive=true)

end # module
