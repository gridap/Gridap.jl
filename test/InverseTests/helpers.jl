using Gridap, Test
using Gridap.FESpaces, Gridap.CellData, Gridap.Arrays

# import Zygote
import Gridap.FESpaces: _compute_cell_ids
import Gridap.CellData: _point_to_cell!, _point_to_cell_cache
import Gridap.Fields: inverse_map
import Gridap.Geometry: get_cell_map
import SparseArrays: sparse


# Degree of freedoms weights matrices computation
function get_dofs_weights_matrix(uh, points)
  Ω = get_triangulation(uh)
  Ug = get_fe_space(uh)
  cell_basis = get_data(get_fe_basis(Ug.space))
  cache = _point_to_cell_cache(KDTreeSearch(), Ω)
  x_to_cell(x) = _point_to_cell!(cache, x)
  point_to_cell = map(x_to_cell, points)
  cell_to_points, _ = make_inverse_table(point_to_cell, num_cells(Ω))
  cell_to_xs = lazy_map(Broadcasting(Reindex(points)), cell_to_points)
  cell_map = get_cell_map(Ω)
  inv_map = lazy_map(inverse_map, cell_map)
  cell_xs_at_ref_space = lazy_map(evaluate, inv_map, cell_to_xs)
  cell_to_weights = lazy_map(evaluate, cell_basis, cell_xs_at_ref_space)
  assemble_dof_weights_matrix(length(points), cell_to_weights, cell_to_points, Ug)
end

# Degree of freedoms weights matrices assembly
function assemble_dof_weights_matrix(npoints, cell_to_weights, cell_to_points, Ug)
  cell_to_dof_ids = get_cell_dof_ids(Ug)
  idxdofs, idxdofspoints, dofentries = Int64[], Int64[], eltype(eltype(cell_to_weights))[]
  idxdirs, idxdirspoints, direntries = Int64[], Int64[], eltype(eltype(cell_to_weights))[]

  function assemble_cell_wisely(weights, ipoints, dof_ids)
    if isempty(weights)
      return
    end

    npts, ndofs = size(weights)
    for ipt in 1:npts
      for idof in 1:ndofs
        if dof_ids[idof] > 0
          push!(idxdofs, dof_ids[idof])
          push!(idxdofspoints, ipoints[ipt])
          push!(dofentries, weights[ipt, idof])
        else
          push!(idxdirs, -dof_ids[idof])
          push!(idxdirspoints, ipoints[ipt])
          push!(direntries, weights[ipt, idof])
        end
      end
    end
  end
  map(assemble_cell_wisely, cell_to_weights, cell_to_points, cell_to_dof_ids)

  ndofs, ndirs = num_free_dofs(Ug), num_dirichlet_dofs(Ug)
  FW = sparse(idxdofspoints, idxdofs, dofentries, npoints, ndofs)
  DW = sparse(idxdirspoints, idxdirs, direntries, npoints, ndirs)
  return FW, DW
end

function interpolate_param_function(func_of_p_and_x, param, Ph)
  func_of_x(x) = func_of_p_and_x(param, x)
  free_values = Vector{eltype(param)}(undef, num_free_dofs(Ph))
  diri_values = Vector{eltype(param)}(undef, num_dirichlet_dofs(Ph))
  interpolate_everywhere!(func_of_x, free_values, diri_values, Ph)
  return free_values
end

function test_gradient(f, x, ϵ=1e-5)
  auto_grad = Zygote.gradient(f, x)[1]
  manu_grad = cdf_gradient(f, x)
  @test maximum(abs, auto_grad - manu_grad) < ϵ
end

function cdf_gradient(f, x, h=1e-4)
  manu_grad = zeros(length(x))
  x⁺, x⁻ = x[:], x[:]
  for i in eachindex(manu_grad)
    x⁺[i] += h
    x⁻[i] -= h

    manu_grad[i] = (f(x⁺) - f(x⁻)) / (2 * h)

    x⁺[i] -= h
    x⁻[i] += h
  end
  manu_grad
end