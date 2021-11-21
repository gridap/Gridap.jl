module RefinedDiscreteModelsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Io

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
@show cell_map = get_cell_map(get_triangulation(model))
@show num_cells = length(cell_map)
@show cell_mask = fill(true, num_cells)
model_ref = newest_vertex_bisection(model, cell_mask)
@test model_ref isa DiscreteModel


end # module
