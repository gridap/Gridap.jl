module FETerms

using Gridap
using Gridap.Helpers

export FETerm
export AffineFETerm
export FESource
export LinearFETerm
export NonLinearFETerm

export setup_cell_jacobian
export setup_cell_residual
export setup_cell_matrix
export setup_cell_vector
export setup_cell_ids

# Interfaces

abstract type FETerm end

"""
Returns an object (e.g. a CellMatrix for single field problems) representing
the contribution to the Jacobian of the given term. Returns nothing if the term
has not contribution to the Jacobian (typically for source terms)
"""
setup_cell_jacobian(t::FETerm,uh,v,du) = @abstractmethod

"""
Returns an object (e.g. a CellVector) representing the contribution to the
residual of the given term. Returns always something.
"""
setup_cell_residual(t::FETerm,uh,v) = @abstractmethod

setup_cell_ids(t::FETerm)::CellNumber = @abstractmethod

abstract type AffineFETerm <: FETerm end

"""
Returns an object (e.g. a CellMatrix) representing the contribution to the
system matrix of the given term. Returns nothing if the term has not
contribution (typically for source terms)
"""
setup_cell_matrix(t::AffineFETerm,v,u) = @abstractmethod

"""
Returns an object (e.g. a CellVector) representing the contribution to the
system rhs of the given term. Returns nothing if the term has not
contribution (typically for linear terms)
"""
setup_cell_vector(t::AffineFETerm,v,uhd) = @abstractmethod

function setup_cell_jacobian(t::AffineFETerm,uh,v,du)
  setup_cell_matrix(t,v,du)
end

# Pretty printing

import Base: show

function show(io::IO,self::FETerm)
  print(io,"$(nameof(typeof(self))) object")
end

function show(io::IO,::MIME"text/plain",self::FETerm)
  show(io,self)
end

# Concrete implementations

struct AffineFETermFromIntegration <: AffineFETerm
  biform::Function
  liform::Function
  trian::Triangulation
  quad::CellQuadrature
end

function AffineFETerm(
  biform::Function,liform::Function,trian::Triangulation,quad::CellQuadrature)
  AffineFETermFromIntegration(biform,liform,trian,quad)
end

function setup_cell_matrix(t::AffineFETermFromIntegration,v,u)
  _v = _restrict_if_needed(v,t.trian)
  _u = _restrict_if_needed(u,t.trian)
  integrate(t.biform(_v,_u),t.trian,t.quad)
end

function setup_cell_vector(t::AffineFETermFromIntegration,v,uhd)
  _v = _restrict_if_needed(v,t.trian)
  _uhd = _restrict_if_needed(uhd,t.trian)
  integrate(t.liform(_v)-t.biform(_v,_uhd),t.trian,t.quad)
end

function setup_cell_residual(t::AffineFETermFromIntegration,uh,v)
  _v = _restrict_if_needed(v,t.trian)
  _uh = _restrict_if_needed(uh,t.trian)
  integrate(t.biform(_v,_uh)-t.liform(_v),t.trian,t.quad)
end

setup_cell_ids(t::AffineFETermFromIntegration) = _setup_cell_ids(t.trian)

struct FESourceFromIntegration <: AffineFETerm
  liform::Function
  trian::Triangulation
  quad::CellQuadrature
end

function FESource(
  liform::Function,
  trian::Triangulation,
  quad::CellQuadrature)
  FESourceFromIntegration(liform,trian,quad)
end

setup_cell_matrix(t::FESourceFromIntegration,v,u) = nothing

function setup_cell_vector(t::FESourceFromIntegration,v,uhd)
  _v = _restrict_if_needed(v,t.trian)
  integrate(t.liform(_v),t.trian,t.quad)
end

function setup_cell_residual(t::FESourceFromIntegration,uh,v)
  _v = _restrict_if_needed(v,t.trian)
  integrate(-t.liform(_v),t.trian,t.quad)
end

setup_cell_ids(t::FESourceFromIntegration) = _setup_cell_ids(t.trian)

struct LinearFETermFromIntegration <: AffineFETerm
  biform::Function
  trian::Triangulation
  quad::CellQuadrature
end

function LinearFETerm(
  biform::Function,trian::Triangulation,quad::CellQuadrature)
  LinearFETermFromIntegration(biform,trian,quad)
end

function setup_cell_matrix(t::LinearFETermFromIntegration,v,u)
  _v = _restrict_if_needed(v,t.trian)
  _u = _restrict_if_needed(u,t.trian)
  integrate(t.biform(_v,_u),t.trian,t.quad)
end

function setup_cell_vector(t::LinearFETermFromIntegration,v,uhd)
  _v = _restrict_if_needed(v,t.trian)
  _uhd = _restrict_if_needed(uhd,t.trian)
  integrate(-t.biform(_v,_uhd),t.trian,t.quad)
end

function setup_cell_residual(t::LinearFETermFromIntegration,uh,v)
  _v = _restrict_if_needed(v,t.trian)
  _uh = _restrict_if_needed(uh,t.trian)
  integrate(t.biform(_v,_uh),t.trian,t.quad)
end

setup_cell_ids(t::LinearFETermFromIntegration) = _setup_cell_ids(t.trian)

struct NonLinearFETerm <: FETerm
  res::Function
  jac::Function
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

# Dealing with several FETerms

function setup_cell_jacobian(uh,v,du,terms::Vararg{<:FETerm})
  m = []
  for term in terms
    r = _jac(term,uh,v,du)
    if r != nothing
      c = _cellids(term)
      _append_matrix_contribution!(m,r,c)
    end
  end
  m
end

function setup_cell_residual(uh,v,terms::Vararg{<:FETerm})
  m = []
  for term in terms
    r = setup_cell_residual(term,uh,v)
    c = _cellids(term)
    _append_vector_contribution!(m,r,c)
  end
  m
end

function setup_cell_matrix(v,u,terms::Vararg{<:AffineFETerm})
  m = []
  for term in terms
    r = _mat(term,v,u)
    if r != nothing
      c = _cellids(term)
      _append_matrix_contribution!(m,r,c)
    end
  end
  m
end

function setup_cell_vector(v,uhd,terms::Vararg{<:AffineFETerm})
  m = []
  for term in terms
    r = _vec(term,v,uhd)
    if r != nothing
      c = _cellids(term)
      _append_vector_contribution!(m,r,c)
    end
  end
  m
end

# Helpers

function _append_vector_contribution!(m,r,c)
  push!(m,(r,c))
end

function _append_vector_contribution!(m,r::SkeletonCellVector,c)
  push!(m,(r.cellvector1,c[1]))
  push!(m,(r.cellvector2,c[2]))
end

function _append_matrix_contribution!(m,r,c)
  push!(m,(r,c,c))
end

function _append_matrix_contribution!(m,r::SkeletonCellMatrix,c)
  push!(m,(r.cellmatrix11,c[1],c[1]))
  push!(m,(r.cellmatrix12,c[1],c[2]))
  push!(m,(r.cellmatrix21,c[2],c[1]))
  push!(m,(r.cellmatrix22,c[2],c[2]))
end

const _jac = setup_cell_jacobian

const _mat = setup_cell_matrix

const _vec = setup_cell_vector

const _cellids = setup_cell_ids

_restrict_if_needed(u,trian) = u

_restrict_if_needed(u,trian::BoundaryTriangulation) = restrict(u,trian)

_restrict_if_needed(u,trian::SkeletonTriangulation) = restrict(u,trian)

_setup_cell_ids(trian) = IdentityCellNumber(Int,ncells(trian))

_setup_cell_ids(trian::BoundaryTriangulation) = trian.descriptor.facet_to_cell

function _setup_cell_ids(trian::SkeletonTriangulation)
  n1 = trian.descriptor1.facet_to_cell
  n2 = trian.descriptor2.facet_to_cell
  (n1, n2)
end

end # module
