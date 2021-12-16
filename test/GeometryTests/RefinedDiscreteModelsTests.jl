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
seed = 5 # Arbitrary
Nsteps = 13
Nsteps_arr = UnitRange(2:Nsteps)
 est = ConstantEst(1.0)
θ = 1.0
write_to_vtk = false
# Uniform refinement
model_refs = build_refined_models(domain, partition, Nsteps, θ, est)
for n = Nsteps
  trian_ref = get_triangulation(model_refs[n])
  if write_to_vtk
    writevtk(trian_ref, "refined$(i)")
  end
  cell_map = get_cell_map(trian_ref)
  node_coords = get_node_coordinates(trian_ref)
  ncoords = length(node_coords)
  # Combinatorial checks for nodes
  if isodd(n)
    a = Integer.(2 * (4^((n-1)/2) + 2^((n-1)/2)) + 1)
  else
    a = Integer.(2^(n/2) + 1)^2
  end
  ncells = length(cell_map)
  @test a == ncoords
  # Combinatorial checks for cells 
  @test ncells == 2^(n + 1)
end
end
