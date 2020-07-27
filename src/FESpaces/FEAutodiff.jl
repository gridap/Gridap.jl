
"""
"""
function autodiff_cell_residual_from_energy(uh_to_cell_energy::Function,uh)
  @assert is_a_fe_function(uh)
  U = get_fe_space(uh)
  cell_u_to_cell_energy = _change_argument_to_cell_u(uh_to_cell_energy,U)
  cell_u = get_cell_values(uh)
  cell_r = autodiff_array_gradient(cell_u_to_cell_energy,cell_u)
  cell_r
end

"""
"""
function autodiff_cell_jacobian_from_energy(uh_to_cell_energy::Function,uh)
  @assert is_a_fe_function(uh)
  U = get_fe_space(uh)
  cell_u_to_cell_energy = _change_argument_to_cell_u(uh_to_cell_energy,U)
  cell_u = get_cell_values(uh)
  cell_j = autodiff_array_hessian(cell_u_to_cell_energy,cell_u)
  cell_j
end

"""
"""
function autodiff_cell_jacobian_from_residual(uh_to_cell_residual::Function,uh)
  @assert is_a_fe_function(uh)
  U = get_fe_space(uh)
  cell_u_to_cell_residual = _change_argument_to_cell_u(uh_to_cell_residual,U)
  cell_u = get_cell_values(uh)
  cell_j = autodiff_array_jacobian(cell_u_to_cell_residual,cell_u)
  cell_j
end

function _change_argument_to_cell_u(uh_to_cell_energy,U)
  function f(cell_u)
    cell_shapefuns = get_cell_shapefuns(U)
    cell_map = get_cell_map(U)
    ref_style = RefStyle(get_cell_dof_basis(U))
    array = lincomb(cell_shapefuns,cell_u)
    uh = GenericCellField(array,cell_map,ref_style)
    cell_e = uh_to_cell_energy(uh)
    cell_e
  end
  f
end

