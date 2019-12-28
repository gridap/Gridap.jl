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

end # module
