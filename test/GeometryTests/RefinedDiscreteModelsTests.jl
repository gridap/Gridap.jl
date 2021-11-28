module RefinedDiscreteModelsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Io

domain = (0,1,0,1)
partition = (1,1)
model = CartesianDiscreteModel(domain,partition)
model = simplexify(model)
cell_map = get_cell_map(get_triangulation(model))
@show num_cells = length(cell_map)
cell_mask = fill(true, num_cells)
model_ref = newest_vertex_bisection(model, cell_mask)
model_ref isa DiscreteModel


end # module
