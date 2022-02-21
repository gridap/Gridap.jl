module NewestVertexBisectionTests

using Test
using Random
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
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

function make_nvb_levels(
  model::DiscreteModel,
  Nsteps::Integer,
  θ::AbstractFloat,
  est::Estimator,
)
  model_refs = Vector{DiscreteModel}(undef, Nsteps)
  cell_map = get_cell_map(get_triangulation(model))
  ncells = length(cell_map)
  η_arr = compute_estimator(est, ncells)
  model_refs[1], buffer = newest_vertex_bisection(model, η_arr; θ = θ)
  buffer = deepcopy(buffer)
  for i = 1:(Nsteps - 1)
    @show i
    cell_map = get_cell_map(get_triangulation(model_refs[i]))
    ncells = length(cell_map)
    η_arr = compute_estimator(est, ncells)
    model_refs[i + 1], buffer =
      newest_vertex_bisection(model_refs[i], buffer, η_arr; θ = θ)
  # Necessary to have each level stored seperately
  buffer = deepcopy(buffer)
  end
  model_refs
end

# For testing only
compute_estimator(est::RandomEst, ncells) = rand(ncells)
compute_estimator(est::ConstantEst, ncells) = fill(est.val, ncells)

domain = (0, 1, 0, 1)
partition = (1, 1) # Initial partition
Nsteps = 12
est = ConstantEst(1.0)
θ = 1.0
uniform_write_to_vtk = false
# Uniform refinement
model = simplexify(CartesianDiscreteModel(domain, partition))
@time model_refs = make_nvb_levels(model, Nsteps, θ, est)
for (n, model_ref) in enumerate(model_refs)
  trian_ref = get_triangulation(model_ref)
  grid = get_grid(model)
  cell_node_ids = get_cell_node_ids(model_ref)
  if uniform_write_to_vtk
    writevtk(trian_ref, "uniform$(string(n, pad=2))")
  end
  cell_map = get_cell_map(trian_ref)
  node_coords = get_node_coordinates(trian_ref)
  ncoords = length(node_coords)
  # Combinatorial checks for nodes
  if isodd(n)
    ncoords_true = Integer.(2 * (4^((n - 1) / 2) + 2^((n - 1) / 2)) + 1)
  else
    ncoords_true = Integer.(2^(n / 2) + 1)^2
  end
  ncells = length(cell_map)
  #@show ncoords
  @test ncoords_true == ncoords
  # Combinatorial checks for cells
  #@test ncells == 2^(n + 1)
end
# Nonuniform refinement. For now only visually checking conformity
#domain = (0, 1, 0, 1)
#partition = (1, 1) # Initial partition
#Nsteps = 20
#seed = 5
#est = RandomEst(seed)
#θ = 0.5
#nonuniform_write_to_vtk = false
#model = simplexify(CartesianDiscreteModel(domain, partition))
#@time model_refs = make_nvb_levels(model, Nsteps, θ, est)
#if nonuniform_write_to_vtk
#  for (n, model_ref) in enumerate(model_refs)
#    trian_ref = get_triangulation(model_ref)
#    writevtk(trian_ref, "nonuniform$(string(n, pad=2))")
#  end
#end

end
