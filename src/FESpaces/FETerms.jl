
"""
    abstract type FETerm <: GridapType end

A `FETerm` is a lazy representation of a summand of a finite element problem. It is not assembled.
See also [FEOperator](@ref).
"""
abstract type FETerm <: GridapType end

"""
Returns an object representing
the contribution to the Jacobian of the given term. Returns nothing if the term
has not contribution to the Jacobian (typically for source terms)
"""
function get_cell_jacobian(t::FETerm,uh,du,v)
  @abstractmethod
end

"""
Returns an object representing the contribution to the
residual of the given term. Returns always something.
"""
function get_cell_residual(t::FETerm,uh,v)
  @abstractmethod
end

"""
"""
function get_cell_id(t::FETerm)
  @abstractmethod
end

"""
"""
function get_cell_jacobian_and_residual(t::FETerm,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  celljac = get_cell_jacobian(t,uh,du,v)
  cellres = get_cell_residual(t,uh,v)
  _setup_jac_and_res(celljac,cellres)
end

function _setup_jac_and_res(celljac,cellres)
  celljacres = pair_arrays(celljac,cellres)
  (celljacres, nothing, nothing)
end

function _setup_jac_and_res(celljac::Nothing,cellres)
  (nothing, nothing, cellres)
end

function _setup_jac_and_res(celljac,cellres::Nothing)
  (nothing, celljac, nothing)
end

function _setup_jac_and_res(celljac::Nothing,cellres::Nothing)
  (nothing, nothing, nothing)
end

"""
See also [FETerm](@ref).
"""
abstract type AffineFETerm <: FETerm end

"""
Returns an object representing the contribution to the
system matrix of the given term. Returns nothing if the term has not
contribution (typically for source terms)
"""
function get_cell_matrix(t::AffineFETerm,u,v)
  @abstractmethod
end

"""
Returns an object (e.g. a CellVector) representing the contribution to the
system rhs of the given term (with Dirichlet bcs included). Returns nothing if the term has not
contribution (typically for linear terms)
"""
function get_cell_vector(t::AffineFETerm,uhd,v)
  @abstractmethod
end

"""
Returns an object (e.g. a CellVector) representing the contribution to the
system rhs of the given term (without Dirichlet bcs). Returns nothing if the term has not
contribution (typically for linear terms)
"""
function get_cell_vector(t::AffineFETerm,v)
  @abstractmethod
end

function get_cell_jacobian(t::AffineFETerm,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  get_cell_matrix(t,du,v)
end

"""
"""
function get_cell_matrix_and_vector(t::AffineFETerm,uhd,u,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  cellmat = get_cell_matrix(t,u,v)
  cellvec = get_cell_vector(t,v)
  cellvals = get_cell_values(uhd,get_cell_id(t))
  _setup_cell_matrix_and_vector(cellmat,cellvec,cellvals)
end

function  _setup_cell_matrix_and_vector(cellmat,cellvec,cellvals)
  cellmatvec = pair_arrays(cellmat,cellvec)
  cellmatvec_with_diri = attach_dirichlet(cellmatvec,cellvals)
  (cellmatvec_with_diri, nothing, nothing)
end

function _setup_cell_matrix_and_vector(cellmat,cellvec::Nothing,cellvals)
  cellmatvec_with_diri = attach_dirichlet(cellmat,cellvals)
  (cellmatvec_with_diri, nothing, nothing)
end

function _setup_cell_matrix_and_vector(cellmat::Nothing,cellvec,cellvals)
  (nothing, nothing, cellvec)
end

function _setup_cell_matrix_and_vector(cellmat::Nothing,cellvec::Nothing,cellvals)
  (nothing, nothing, nothing)
end

"""
See also [FETerm](@ref).
"""
abstract type LinearFETerm <: AffineFETerm end

"""
"""
abstract type FESource <: AffineFETerm end

function get_cell_matrix(t::FESource,u,v)
  nothing
end

# Working with a collection of FETerms

"""
"""
function collect_cell_jacobian(uh,du,v,terms)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  w = []
  r = []
  c = []
  for term in terms
    cellvals = get_cell_jacobian(term,uh,du,v)
    cellids = get_cell_id(term)
    _push_matrix_contribution!(w,r,c,cellvals,cellids)
  end
  (w,r,c)
end

"""
"""
function collect_cell_matrix(u,v,terms)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  w = []
  r = []
  c = []
  for term in terms
    cellvals = get_cell_matrix(term,u,v)
    cellids = get_cell_id(term)
    _push_matrix_contribution!(w,r,c,cellvals,cellids)
  end
  (w,r,c)
end

"""
"""
function collect_cell_residual(uh,v,terms)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  w = []
  r = []
  for term in terms
    cellvals = get_cell_residual(term,uh,v)
    cellids = get_cell_id(term)
    _push_vector_contribution!(w,r,cellvals,cellids)
  end
  (w,r)
end

"""
"""
function collect_cell_vector(uhd,v,terms)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_function(uhd)
  w = []
  r = []
  for term in terms
    cellvals = get_cell_vector(term,uhd,v)
    cellids = get_cell_id(term)
    _push_vector_contribution!(w,r,cellvals,cellids)
  end
  (w,r)
end

"""
"""
function collect_cell_matrix_and_vector(uhd,u,v,terms)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  @assert is_a_fe_function(uhd)

  matvecdata = ([],[],[])
  matdata = ([],[],[])
  vecdata = ([],[])

  for term in terms
    cellmatvec, cellmat, cellvec = get_cell_matrix_and_vector(term,uhd,u,v)
    cellids = get_cell_id(term)
    _push_matrix_contribution!(matvecdata...,cellmatvec,cellids)
    _push_matrix_contribution!(matdata...,cellmat,cellids)
    _push_vector_contribution!(vecdata...,cellvec,cellids)
  end

  (matvecdata, matdata, vecdata)
end

"""
"""
function collect_cell_jacobian_and_residual(uh,du,v,terms)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)

  matvecdata = ([],[],[])
  matdata = ([],[],[])
  vecdata = ([],[])

  for term in terms
    cellmatvec, cellmat, cellvec = get_cell_jacobian_and_residual(term,uh,du,v)
    cellids = get_cell_id(term)
    _push_matrix_contribution!(matvecdata...,cellmatvec,cellids)
    _push_matrix_contribution!(matdata...,cellmat,cellids)
    _push_vector_contribution!(vecdata...,cellvec,cellids)
  end

  (matvecdata, matdata, vecdata)
end

function _push_matrix_contribution!(w,r,c,cellvals,cellids)
  push!(w,cellvals)
  push!(r,cellids)
  push!(c,cellids)
  nothing
end

function _push_matrix_contribution!(w,r,c,cellvals::Nothing,cellids)
  nothing
end

function _push_vector_contribution!(v,r,cellvals,cellids)
  push!(v,cellvals)
  push!(r,cellids)
  nothing
end

function _push_vector_contribution!(v,r,cellvals::Nothing,cellids)
  nothing
end

# Concrete implementations

struct AffineFETermFromIntegration <: AffineFETerm
  biform::Function
  liform::Function
  trian::Triangulation
  quad::CellQuadrature
end

"""
"""
function AffineFETerm(
  biform::Function,liform::Function,trian::Triangulation,quad::CellQuadrature)
  AffineFETermFromIntegration(biform,liform,trian,quad)
end

function get_cell_matrix(t::AffineFETermFromIntegration,u,v)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  integrate(t.biform(u,v),t.quad)
end

function get_cell_vector(t::AffineFETermFromIntegration,uhd,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  integrate(t.liform(v)-t.biform(uhd,v),t.quad)
end

function get_cell_vector(t::AffineFETermFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  integrate(t.liform(v),t.quad)
end

function get_cell_residual(t::AffineFETermFromIntegration,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  integrate(t.biform(uh,v)-t.liform(v),t.quad)
end

function get_cell_id(t::AffineFETermFromIntegration)
  get_cell_id(t.trian)
end

struct FESourceFromIntegration <: FESource
  liform::Function
  trian::Triangulation
  quad::CellQuadrature
end

"""
"""
function FESource(
  liform::Function,
  trian::Triangulation,
  quad::CellQuadrature)
  FESourceFromIntegration(liform,trian,quad)
end

function get_cell_matrix(t::FESourceFromIntegration,u,v)
  nothing
end

function get_cell_vector(t::FESourceFromIntegration,uhd,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  integrate(t.liform(v),t.quad)
end

function get_cell_vector(t::FESourceFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  integrate(t.liform(v),t.quad)
end

function get_cell_residual(t::FESourceFromIntegration,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  integrate(-t.liform(v),t.quad)
end

function get_cell_id(t::FESourceFromIntegration)
  get_cell_id(t.trian)
end

struct LinearFETermFromIntegration <: LinearFETerm
  biform::Function
  trian::Triangulation
  quad::CellQuadrature
end

"""
"""
function LinearFETerm(
  biform::Function,trian::Triangulation,quad::CellQuadrature)
  LinearFETermFromIntegration(biform,trian,quad)
end

function get_cell_matrix(t::LinearFETermFromIntegration,u,v)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  integrate(t.biform(u,v),t.quad)
end

function get_cell_vector(t::LinearFETermFromIntegration,uhd,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  integrate(-t.biform(uhd,v),t.quad)
end

function get_cell_vector(t::LinearFETermFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  nothing
end

function get_cell_residual(t::LinearFETermFromIntegration,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  integrate(t.biform(uh,v),t.quad)
end

function get_cell_id(t::LinearFETermFromIntegration)
  get_cell_id(t.trian)
end

struct NonlinearFETerm <: FETerm
  res::Function
  jac::Function
  trian::Triangulation
  quad::CellQuadrature
end

"""
"""
function FETerm(
  res::Function, jac::Function, trian::Triangulation, quad::CellQuadrature)
  NonlinearFETerm(res,jac,trian,quad)
end

function get_cell_jacobian(t::NonlinearFETerm,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  integrate(t.jac(uh,du,v),t.quad)
end

function get_cell_residual(t::NonlinearFETerm,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  integrate(t.res(uh,v),t.quad)
end

function get_cell_id(t::NonlinearFETerm)
  get_cell_id(t.trian)
end

"""
"""
struct AffineFETermFromCellMatVec <: AffineFETerm
  matvecfun::Function
  trian::Triangulation
end

function get_cell_matrix(t::AffineFETermFromCellMatVec,u,v)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  cellmatvec = t.matvecfun(u,v)
  cellmat, _ = unpair_arrays(cellmatvec)
  cellmat
end

function get_cell_vector(t::AffineFETermFromCellMatVec,uhd,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  trial = TrialFESpace(get_fe_space(uhd))
  u = get_cell_basis(trial)
  _cellvals = get_cell_values(uhd,get_cell_id(t.trian))
  cellmatvec = t.matvecfun(u,v)
  cellmatvec_with_diri = attach_dirichlet(cellmatvec,_cellvals)
  _, cellvec_with_diri = unpair_arrays(cellmatvec_with_diri)
  cellvec_with_diri
end

function get_cell_vector(t::AffineFETermFromCellMatVec,v)
  @assert is_a_fe_cell_basis(v)
  @unreachable "This function cannot be implemented for $(typeof(t)) objects."
end

function get_cell_id(t::AffineFETermFromCellMatVec)
  get_cell_id(t.trian)
end

function get_cell_matrix_and_vector(t::AffineFETermFromCellMatVec,uhd,u,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  _cellvals = get_cell_values(uhd,get_cell_id(t.trian))
  cellmatvec = t.matvecfun(u,v)
  cellmatvec_with_diri = attach_dirichlet(cellmatvec,_cellvals)
  (cellmatvec_with_diri, nothing, nothing)
end

"""
"""
struct FETermFromCellJacRes <: FETerm
  jacresfun::Function
  trian::Triangulation
end

function get_cell_jacobian(t::FETermFromCellJacRes,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  celljacres = t.jacresfun(uh,du,v)
  celljac, _ = unpair_arrays(celljacres)
  celljac
end

function get_cell_residual(t::FETermFromCellJacRes,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  trial = TrialFESpace(get_fe_space(uh))
  du = get_cell_basis(trial)
  celljacres = t.jacresfun(uh,du,v)
  _, cellres = unpair_arrays(celljacres)
  cellres
end

function get_cell_id(t::FETermFromCellJacRes)
  get_cell_id(t.trian)
end

function get_cell_jacobian_and_residual(t::FETermFromCellJacRes,uh,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  celljacres = t.jacresfun(uh,du,v)
  (celljacres, nothing, nothing)
end
