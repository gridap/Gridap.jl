

function gradient(f::Function,uh::FEFunction)
  fuh = f(uh)
  _gradient(f,uh,fuh)
end

function _gradient(f,uh,fuh::AbstractArray)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.
  
  Make sure that you are using a LebesgueMeasure instead of a CellQuadrature to perform integration.
  """
end

function _gradient(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = get_cell_id(trian)
    cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function jacobian(f::Function,uh::FEFunction)
  fuh = f(uh)
  _jacobian(f,uh,fuh)
end

function _jacobian(f,uh,fuh::AbstractArray)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.
  
  Make sure that you are using a LebesgueMeasure instead of a CellQuadrature to perform integration.
  """
end

function _jacobian(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = get_cell_id(trian)
    cell_grad = autodiff_array_jacobian(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function hessian(f::Function,uh::FEFunction)
  fuh = f(uh)
  _hessian(f,uh,fuh)
end

function _hessian(f,uh,fuh::AbstractArray)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.
  
  Make sure that you are using a LebesgueMeasure instead of a CellQuadrature to perform integration.
  """
end

function _hessian(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = get_cell_id(trian)
    cell_grad = autodiff_array_hessian(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function _change_argument(f,trian,uh)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = CellField(U,cell_u)
    cell_grad = f(cf)
    get_contribution(cell_grad,trian)
  end
  g
end


