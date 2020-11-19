
"""
"""
function autodiff_cell_residual_from_energy(
  uh_to_cell_energy::Function,
  uh::FEFunction,
  cell_ids=_default_cell_ids(uh))

  U = get_fe_space(uh)
  cell_u_to_cell_energy = _change_argument_to_cell_u(uh_to_cell_energy,U)
  cell_u = get_cell_dof_values(uh)
  cell_r = autodiff_array_gradient(cell_u_to_cell_energy,cell_u,cell_ids)
  cell_r
end

"""
"""
function autodiff_cell_jacobian_from_energy(
  uh_to_cell_energy::Function,
  uh::FEFunction,
  cell_ids=_default_cell_ids(uh))

  U = get_fe_space(uh)
  cell_u_to_cell_energy = _change_argument_to_cell_u(uh_to_cell_energy,U)
  cell_u = get_cell_dof_values(uh)
  cell_j = autodiff_array_hessian(cell_u_to_cell_energy,cell_u,cell_ids)
  cell_j
end

"""
"""
function autodiff_cell_jacobian_from_residual(
  uh_to_cell_residual::Function,
  uh::FEFunction,
  cell_ids=_default_cell_ids(uh))

  U = get_fe_space(uh)
  cell_u_to_cell_residual = _change_argument_to_cell_u(uh_to_cell_residual,U)
  cell_u = get_cell_dof_values(uh)
  cell_j = autodiff_array_jacobian(cell_u_to_cell_residual,cell_u,cell_ids)
  cell_j
end

function _default_cell_ids(uh)
  cell_u = get_cell_dof_values(uh)
  cell_ids = IdentityVector(length(cell_u))
end

function _change_argument_to_cell_u(uh_to_cell_energy,U)
  function f(cell_u)
    uh = CellField(U,cell_u)
    cell_e = uh_to_cell_energy(uh)
    get_array(cell_e)
  end
  f
end

