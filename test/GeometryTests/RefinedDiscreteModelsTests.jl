module RefinedDiscreteModelsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Visualization

domain = (0,1,0,1)
partition = (1,1)
model = CartesianDiscreteModel(domain,partition)
model_ref = simplexify(model)
cell_map = get_cell_map(get_triangulation(model))
num_cells = length(cell_map)
cell_mask = fill(true, num_cells)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, true)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
writevtk(Triangulation(model_ref), "refined6")
#for i in 1:3
#    @show i
#    @show model_ref
#    #writevtk(Triangulation(model_ref), "refined$(i)")
#    model_ref = newest_vertex_bisection(model_ref, cell_mask, false)
#end
#model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
#model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
#@show model_ref =  newest_vertex_bisection(model_ref, cell_mask, true)
#@show model_ref =  newest_vertex_bisection(model_ref, cell_mask, true)
@test model_ref isa DiscreteModel
end # module
