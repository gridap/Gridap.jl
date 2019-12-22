module UnstructuredDiscreteModelsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock

m = DiscreteModelMock()
g = get_grid(m)

model = UnstructuredDiscreteModel(g)
test_discrete_model(model)

end # module
