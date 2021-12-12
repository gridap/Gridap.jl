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
η_arr = rand(num_cells)
Nsteps = 15
model_refs = Vector{DiscreteModel}(undef, Nsteps)
sort_flag = true
θ = 0.5
model_refs[1] = newest_vertex_bisection(model, η_arr; sort_flag=sort_flag, θ=θ)
for i in 1:Nsteps-1
    cell_map = get_cell_map(get_triangulation(model_refs[i]))
    num_cells = length(cell_map)
    η_arr = rand(num_cells)
    @show i
    writevtk(Triangulation(model_refs[i]), "refined$(i)")
    model_refs[i + 1] = newest_vertex_bisection(model_refs[i], η_arr, sort_flag=sort_flag, θ=θ)
end

#model_ref =  newest_vertex_bisection(model_ref, η_arr, false)
#model_ref =  newest_vertex_bisection(model_ref, η_arr, false)
#@show model_ref =  newest_vertex_bisection(model_ref, η_arr, true)
#@show model_ref =  newest_vertex_bisection(model_ref, η_arr, true)
@test model_refs[end] isa DiscreteModel
end # module
