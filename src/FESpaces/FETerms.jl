module FETerms

using Gridap
using Gridap.Helpers

abstract type FETerm end

setup_cell_jacobian(t::FETerm,uh,v,du)::CellMatrix = @abstractmethod

setup_cell_residual(t::FETerm,uh,v)::CellVector = @abstractmethod

setup_cell_ids(t::FETerm)::CellNumber = @abstractmethod

abstract type FESource <: FETerm end

setup_cell_vector(t::FESource,v,uhd)::CellVector = @abstractmethod

setup_cell_jacobian(t::FESource,uh,v,du) = nothing

setup_cell_residual(t::FESource,uh,v) = - setup_cell_vector(t,v)

struct LinearFETerm <: FETerm
  biform::Function
  liform::Function
  trian::Triangulation
  quad::CellQuadrature
end

function setup_cell_matrix(t::LinearFETerm,v,u)
  _v = _restrict_if_needed(v,t.trian)
  _u = _restrict_if_needed(u,t.trian)
  integrate(t.biform(_v,_u),t.trian,t.quad)
end

function setup_cell_vector(t::LinearFETerm,v,uhd)
  _v = _restrict_if_needed(v,t.trian)
  _uhd = _restrict_if_needed(uhd,t.trian)
  integrate(t.liform(_v)-t.biform(_v,_uhd),t.trian,t.quad)
end

function setup_cell_jacobian(t::LinearFETerm,uh,v,du)
  setup_cell_matrix(t,v,du)
end

function setup_cell_residual(t::LinearFETerm,uh,v)
  _v = _restrict_if_needed(v,t.trian)
  _uh = _restrict_if_needed(uh,t.trian)
  integrate(t.biform(_v,_uh)-t.liform(_v),t.trian,t.quad)
end

setup_cell_ids(t::LinearFETerm) = _setup_cell_ids(t.trian)

struct NonLinearFETerm <: FETerm
  jac::Function
  res::Function
  trian::Triangulation
  quad::CellQuadrature
end

function setup_cell_jacobian(t::NonLinearFETerm,uh,v,du)
  _v = _restrict_if_needed(v,t.trian)
  _uh = _restrict_if_needed(uh,t.trian)
  _du = _restrict_if_needed(du,t.trian)
  integrate(t.jac(_uh,_v,_du),t.trian,t.quad)
end

function setup_cell_residual(t::NonLinearFETerm,uh,v)
  _v = _restrict_if_needed(v,t.trian)
  _uh = _restrict_if_needed(uh,t.trian)
  integrate(t.res(_uh,_v),t.trian,t.quad)
end

setup_cell_ids(t::NonLinearFETerm) = _setup_cell_ids(t.trian)

struct WeakFESource <: FESource
  liform::Function
  trian::Triangulation
  quad::CellQuadrature
end

function FESource(
  liform::Function,trian::Triangulation,quad::CellQuadrature)
  WeakFESource(liform,trian,quad)
end

function setup_cell_vector(t::WeakFESource,v,uhd)
  _v = _restrict_if_needed(v,t.trian)
  integrate(t.liform(_v),t.trian,t.quad)
end

setup_cell_ids(t::WeakFESource) = _setup_cell_ids(t.trian)

_restrict_if_needed(u,trian) = u

_restrict_if_needed(u,trian::BoundaryTriangulation) = restrict(u,trian)

_setup_cell_ids(trian) = IdentityCellNumber(Int,ncells(trian))

_setup_cell_ids(trian::BoundaryTriangulation) = trian.descriptor.facet_to_cell

end # module
