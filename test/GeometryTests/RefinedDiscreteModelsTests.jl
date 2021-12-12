module RefinedDiscreteModelsTests

using Test
using Random
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Visualization


# For testing only
abstract type Estimator end

struct ConstantEst <: Estimator
  val::Float64
end

struct RandomEst <: Estimator
  function RandomEst(seed)
    Random.seed!(seed)
    new()
  end
end

# For testing only
compute_estimator(est::RandomEst, ncells) = rand(ncells)
compute_estimator(est::ConstantEst, ncells) = fill(est.val, ncells)



domain = (0, 1, 0, 1)
partition = (1, 1) # Initial partition
seed = 5 # Arbitrary
model = CartesianDiscreteModel(domain, partition)
model = simplexify(model)
cell_map = get_cell_map(get_triangulation(model))
ncells = length(cell_map)
est = RandomEst(seed)
#est = ConstantEst(1.0)
η_arr = compute_estimator(est, ncells)
Nsteps = 14
model_refs = Vector{DiscreteModel}(undef, Nsteps)
sort_flag = true
θ = 0.6
model_refs[1] = newest_vertex_bisection(model, η_arr; sort_flag = sort_flag, θ = θ)
sort_flag = false
for i = 1:(Nsteps - 1)
  cell_map = get_cell_map(get_triangulation(model_refs[i]))
  ncells = length(cell_map)
  η_arr = compute_estimator(est, ncells)
  #writevtk(Triangulation(model_refs[i]), "refined$(i)")
  model_refs[i + 1] =
    newest_vertex_bisection(model_refs[i], η_arr, sort_flag = sort_flag, θ = θ)
end

@test model_refs[end] isa DiscreteModel
end # module
