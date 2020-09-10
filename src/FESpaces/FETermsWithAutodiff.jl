"""
"""
function FETerm(
  res::Function, trian::Triangulation, quad::CellQuadrature)
  NonlinearFETermWithAutodiff(res,trian,quad)
end

struct NonlinearFETermWithAutodiff <: FETerm
  res::Function
  trian::Triangulation
  quad::CellQuadrature
end

function get_cell_residual(t::NonlinearFETermWithAutodiff,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  integrate(t.res(uh,v),t.quad)
end

function get_cell_jacobian(t::NonlinearFETermWithAutodiff,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  function uh_to_cell_residual(uh)
    integrate(t.res(uh,v),t.quad)
  end
  cell_j = autodiff_cell_jacobian_from_residual(uh_to_cell_residual,uh,get_cell_id(t))
  cell_j
end

function get_cell_id(t::NonlinearFETermWithAutodiff)
  get_cell_id(t.trian)
end

struct FEEnergy <: FETerm
  ener::Function
  trian::Triangulation
  quad::CellQuadrature
end

function get_cell_residual(t::FEEnergy,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  function uh_to_cell_energy(uh)
    integrate(t.ener(uh),t.quad)
  end
  cell_r = autodiff_cell_residual_from_energy(uh_to_cell_energy,uh,get_cell_id(t))
  cell_r
end

function get_cell_jacobian(t::FEEnergy,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  function uh_to_cell_energy(uh)
    integrate(t.ener(uh),t.quad)
  end
  cell_j = autodiff_cell_jacobian_from_energy(uh_to_cell_energy,uh,get_cell_id(t))
  cell_j
end

function get_cell_id(t::FEEnergy)
  get_cell_id(t.trian)
end
