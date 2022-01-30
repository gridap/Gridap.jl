module RefinedDiscreteModelsTests

using Test
using Random
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: ConstantEst, RandomEst
using Gridap.Geometry: build_refined_models
using Gridap.Visualization
#using TimerOutputs

domain = (0, 1, 0, 1)
partition = (1, 1) # Initial partition
Nsteps = 12
est = ConstantEst(1.0)
θ = 1.0
uniform_write_to_vtk = false
# Uniform refinement
model = simplexify(CartesianDiscreteModel(domain, partition))
model_refs = build_refined_models(model, Nsteps, θ, est)
for (n, model_ref) in enumerate(model_refs)
  trian_ref = get_triangulation(model_ref)
  if uniform_write_to_vtk
    writevtk(trian_ref, "uniform$(n)")
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
  @test ncoords_true == ncoords
  # Combinatorial checks for cells
  @test ncells == 2^(n + 1)
end
# Nonuniform refinement. For now only visually checking conformity
#domain = (0, 1, 0, 1)
#partition = (1, 1) # Initial partition
#Nsteps = 13
#seed = 5
#est = RandomEst(seed)
#θ = 0.5
#nonuniform_write_to_vtk = true
#model = simplexify(CartesianDiscreteModel(domain, partition))
#model_refs = build_refined_models(model, Nsteps, θ, est)
#if nonuniform_write_to_vtk
#  for (n, model_ref) in enumerate(model_refs)
#    trian_ref = get_triangulation(model_ref)
#    writevtk(trian_ref, "nonuniform$(n)")
#  end
#end

end
