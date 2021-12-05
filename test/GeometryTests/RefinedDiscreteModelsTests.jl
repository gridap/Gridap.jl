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
model = simplexify(model)
cell_map = get_cell_map(get_triangulation(model))
num_cells = length(cell_map)
cell_mask = fill(true, num_cells)
Nsteps = 30
model_refs = Vector{DiscreteModel}(undef, Nsteps)
model_refs[1] = newest_vertex_bisection(model, cell_mask, true)
for i in 1:Nsteps-1
    @show i
    model_refs[i + 1] = newest_vertex_bisection(model_refs[i], cell_mask, false)
end
writevtk(Triangulation(model_refs[end]), "refined")

#model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
#model_ref =  newest_vertex_bisection(model_ref, cell_mask, false)
#@show model_ref =  newest_vertex_bisection(model_ref, cell_mask, true)
#@show model_ref =  newest_vertex_bisection(model_ref, cell_mask, true)
@test model_refs[end] isa DiscreteModel
end # module
