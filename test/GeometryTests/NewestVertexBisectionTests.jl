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
    cell_map = get_cell_map(get_triangulation(model_refs[i]))
    ncells = length(cell_map)
    η_arr = compute_estimator(est, ncells)
    model_refs[i + 1], buffer =
      newest_vertex_bisection(model_refs[i], buffer, η_arr; θ = θ)
  # Necessary to have each level stored separately
  buffer = deepcopy(buffer)
  end
  model_refs
end

function d_get_num_boundary_labels(labels, d)
  entities = labels.d_to_dface_to_entity[d + 1]
  names = labels.tag_to_name[entities]
  bdry_names = filter(name -> name != "interior", names)
  length(bdry_names)
end


# For testing only
compute_estimator(est::RandomEst, ncells) = rand(ncells)
compute_estimator(est::ConstantEst, ncells) = fill(est.val, ncells)

domain = (0, 1, 0, 1)
partition = (1, 1) # Initial partition
Nsteps = 10
est = ConstantEst(1.0)
θ = 1.0
uniform_write_to_vtk = false
# Uniform refinement
let model = simplexify(CartesianDiscreteModel(domain, partition))
  model_refs = make_nvb_levels(model, Nsteps, θ, est)
  for (n, model_ref) in enumerate(model_refs)
    trian_ref = get_triangulation(model_ref)
    if uniform_write_to_vtk
      writevtk(trian_ref, "uniform$(string(n, pad=2))")
    end
    # Combinatorial checks for cells
    let cell_map = get_cell_map(trian_ref)
      ncells = length(cell_map)
      @test ncells == 2^(n + 1)
    end
    # Combinatorial checks for nodes
    let node_coords = get_node_coordinates(trian_ref)
      ncoords = length(node_coords)
      if isodd(n)
        ncoords_true = Int(2 * (4^((n - 1) / 2) + 2^((n - 1) / 2)) + 1)
      else
        ncoords_true = Int(2^(n / 2) + 1)^2
      end
      @test ncoords_true == ncoords
    end
    # Test labels are propagating properly
    let labels = get_face_labeling(model_ref)
      let d = 0
        @test Int(4*2^(floor(n / 2))) == d_get_num_boundary_labels(labels, d)
      end
      # Combinatorial check for face_labels (edge)
      let d = 1
        @test Int(4*2^(floor(n / 2))) == d_get_num_boundary_labels(labels, d)
      end
      # Combinatorial check for cell labels
      let d = 2
        @test 0 == d_get_num_boundary_labels(labels, d)
      end
    end
  end
end

# Nonuniform refinement. For now only visually checking conformity
#domain = (0, 1, 0, 1)
#partition = (1, 1) # Initial partition
#Nsteps = 15
#seed = 5
#est = RandomEst(seed)
#θ = 0.5
#nonuniform_write_to_vtk = true
#model = simplexify(CartesianDiscreteModel(domain, partition))
#@time model_refs = make_nvb_levels(model, Nsteps, θ, est)
#if nonuniform_write_to_vtk
#  for (n, model_ref) in enumerate(model_refs)
#    trian_ref = get_triangulation(model_ref)
#    writevtk(trian_ref, "nonuniform$(string(n, pad=2))")
#  end
#end

end
