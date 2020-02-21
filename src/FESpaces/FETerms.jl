
"""
"""
abstract type FETerm <: GridapType end

"""
Returns an object representing
the contribution to the Jacobian of the given term. Returns nothing if the term
has not contribution to the Jacobian (typically for source terms)
"""
function get_cell_jacobian(t::FETerm,uh,v,du)
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
function get_cell_jacobian_and_residual(t::FETerm,uh,v,du)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  celljac = get_cell_jacobian(t,uh,v,du)
  cellres = get_cell_residual(t,uh,v)
  _setup_jac_and_res(celljac,cellres)
end

function _setup_jac_and_res(celljac,cellres)
  celljacres = pair_arrays(celljac,cellres)
  (celljacres, nothing, nothing)
end

function _setup_jac_and_res(celljac::SkeletonCellMatrix,cellres::SkeletonCellVector)
  #TODO for the moment do not pair on the skeleton
  (nothing, celljac, cellres)
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
"""
abstract type AffineFETerm <: FETerm end

"""
Returns an object representing the contribution to the
system matrix of the given term. Returns nothing if the term has not
contribution (typically for source terms)
"""
function get_cell_matrix(t::AffineFETerm,v,u)
  @abstractmethod
end

"""
Returns an object (e.g. a CellVector) representing the contribution to the
system rhs of the given term. Returns nothing if the term has not
contribution (typically for linear terms)
"""
function get_cell_vector(t::AffineFETerm,v,uhd)
  @abstractmethod
end

function get_cell_vector(t::AffineFETerm,v)
  @abstractmethod
end

"""
"""
function get_cell_values(t::FETerm,uhd)
  @abstractmethod
end

function get_cell_jacobian(t::AffineFETerm,uh,v,du)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  get_cell_matrix(t,v,du)
end

"""
"""
function get_cell_matrix_and_vector(t::AffineFETerm,v,u,uhd)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  cellmat = get_cell_matrix(t,v,u)
  cellvec = _get_cell_vector_tmp_hack(cellmat,t,v,uhd) #TODO
  cellvals = get_cell_values(t,uhd)
  _setup_cell_matrix_and_vector(cellmat,cellvec,cellvals)
end

function _get_cell_vector_tmp_hack(cellmat,t,v,uhd) #TODO
  cellvec = get_cell_vector(t,v)
  cellvec
end

function _get_cell_vector_tmp_hack(cellmat::SkeletonCellMatrix,t,v,uhd) #TODO
  cellvec = get_cell_vector(t,v,uhd)
  cellvec
end

function  _setup_cell_matrix_and_vector(cellmat,cellvec,cellvals)
  cellmatvec = pair_arrays(cellmat,cellvec)
  cellmatvec_with_diri = attach_dirichlet_bcs(cellmatvec,cellvals)
  (cellmatvec_with_diri, nothing, nothing)
end

function  _setup_cell_matrix_and_vector(cellmat::SkeletonCellMatrix,cellvec::SkeletonCellVector,cellvals)
  # TODO for the moment do not pair quantities on the skeleton and assume that cellvec has dirichlet bcs
  (nothing, cellmat, cellvec)
end

function _setup_cell_matrix_and_vector(cellmat,cellvec::Nothing,cellvals)
  cellmatvec_with_diri = attach_dirichlet_bcs(cellmat,cellvals)
  (cellmatvec_with_diri, nothing, nothing)
end

function _setup_cell_matrix_and_vector(cellmat::Nothing,cellvec,cellvals)
  (nothing, nothing, cellvec)
end

function _setup_cell_matrix_and_vector(cellmat::Nothing,cellvec::Nothing,cellvals)
  (nothing, nothing, nothing)
end

"""
"""
abstract type LinearFETerm <: AffineFETerm end

"""
"""
abstract type FESource <: AffineFETerm end

function get_cell_matrix(t::FESource,v,u)
  nothing
end

# Working with a collection of FETerms

"""
"""
function collect_cell_jacobian(uh,v,du,terms)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  w = []
  r = []
  c = []
  for term in terms
    cellvals = get_cell_jacobian(term,uh,v,du)
    cellids = get_cell_id(term)
    _push_matrix_contribution!(w,r,c,cellvals,cellids)
  end
  (w,r,c)
end

"""
"""
function collect_cell_matrix(v,u,terms)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  w = []
  r = []
  c = []
  for term in terms
    cellvals = get_cell_matrix(term,v,u)
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
function collect_cell_vector(v,uhd,terms)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_function(uhd)
  w = []
  r = []
  for term in terms
    cellvals = get_cell_vector(term,v,uhd)
    cellids = get_cell_id(term)
    _push_vector_contribution!(w,r,cellvals,cellids)
  end
  (w,r)
end

"""
"""
function collect_cell_matrix_and_vector(v,u,uhd,terms)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  @assert is_a_fe_function(uhd)

  matvecdata = ([],[],[])
  matdata = ([],[],[])
  vecdata = ([],[])

  for term in terms
    cellmatvec, cellmat, cellvec = get_cell_matrix_and_vector(term,v,u,uhd)
    cellids = get_cell_id(term)
    _push_matrix_contribution!(matvecdata...,cellmatvec,cellids)
    _push_matrix_contribution!(matdata...,cellmat,cellids)
    _push_vector_contribution!(vecdata...,cellvec,cellids)
  end

  (matvecdata, matdata, vecdata)
end

"""
"""
function collect_cell_jacobian_and_residual(uh,v,du,terms)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)

  matvecdata = ([],[],[])
  matdata = ([],[],[])
  vecdata = ([],[])

  for term in terms
    cellmatvec, cellmat, cellvec = get_cell_jacobian_and_residual(term,uh,v,du)
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

function _push_matrix_contribution!(w,r,c,cellvals::SkeletonCellMatrix,cellids::SkeletonPair)
  push!(w,cellvals.ll)
  push!(w,cellvals.lr)
  push!(w,cellvals.rl)
  push!(w,cellvals.rr)
  push!(r,cellids.left); push!(c,cellids.left)
  push!(r,cellids.left); push!(c,cellids.right)
  push!(r,cellids.right); push!(c,cellids.left)
  push!(r,cellids.right); push!(c,cellids.right)
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

function _push_vector_contribution!(v,r,cellvals::SkeletonCellVector,cellids::SkeletonPair)
  push!(v,cellvals.left)
  push!(v,cellvals.right)
  push!(r,cellids.left)
  push!(r,cellids.right)
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

function get_cell_matrix(t::AffineFETermFromIntegration,v,u)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  _v = restrict(v,t.trian)
  _u = restrict(u,t.trian)
  integrate(t.biform(_v,_u),t.trian,t.quad)
end

function get_cell_vector(t::AffineFETermFromIntegration,v,uhd)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  _uhd = restrict(uhd,t.trian)
  integrate(t.liform(_v)-t.biform(_v,_uhd),t.trian,t.quad)
end

function get_cell_vector(t::AffineFETermFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  integrate(t.liform(_v),t.trian,t.quad)
end

function get_cell_residual(t::AffineFETermFromIntegration,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  _uh = restrict(uh,t.trian)
  integrate(t.biform(_v,_uh)-t.liform(_v),t.trian,t.quad)
end

function get_cell_id(t::AffineFETermFromIntegration)
  get_cell_id(t.trian)
end

function get_cell_values(t::AffineFETermFromIntegration,uhd)
  @assert is_a_fe_function(uhd)
  cellvals = get_cell_values(uhd)
  reindex(cellvals,t.trian)
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

function get_cell_matrix(t::FESourceFromIntegration,v,u)
  nothing
end

function get_cell_vector(t::FESourceFromIntegration,v,uhd)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  integrate(t.liform(_v),t.trian,t.quad)
end

function get_cell_vector(t::FESourceFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  integrate(t.liform(_v),t.trian,t.quad)
end

function get_cell_residual(t::FESourceFromIntegration,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  integrate(-t.liform(_v),t.trian,t.quad)
end

function get_cell_id(t::FESourceFromIntegration)
  get_cell_id(t.trian)
end

function get_cell_values(t::FESourceFromIntegration,uhd)
  @assert is_a_fe_function(uhd)
  cellvals = get_cell_values(uhd)
  reindex(cellvals,t.trian)
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

function get_cell_matrix(t::LinearFETermFromIntegration,v,u)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  _v = restrict(v,t.trian)
  _u = restrict(u,t.trian)
  integrate(t.biform(_v,_u),t.trian,t.quad)
end

function get_cell_vector(t::LinearFETermFromIntegration,v,uhd)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  _uhd = restrict(uhd,t.trian)
  integrate(-t.biform(_v,_uhd),t.trian,t.quad)
end

function get_cell_vector(t::LinearFETermFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  nothing
end

function get_cell_residual(t::LinearFETermFromIntegration,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  _uh = restrict(uh,t.trian)
  integrate(t.biform(_v,_uh),t.trian,t.quad)
end

function get_cell_id(t::LinearFETermFromIntegration) 
  get_cell_id(t.trian)
end

function get_cell_values(t::LinearFETermFromIntegration,uhd)
  @assert is_a_fe_function(uhd)
  cellvals = get_cell_values(uhd)
  reindex(cellvals,t.trian)
end

struct NonLinearFETerm <: FETerm
  res::Function
  jac::Function
  trian::Triangulation
  quad::CellQuadrature
end

"""
"""
function FETerm(
  res::Function, jac::Function, trian::Triangulation, quad::CellQuadrature)
  NonLinearFETerm(res,jac,trian,quad)
end

function get_cell_jacobian(t::NonLinearFETerm,uh,v,du)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  _v = restrict(v,t.trian)
  _uh = restrict(uh,t.trian)
  _du = restrict(du,t.trian)
  integrate(t.jac(_uh,_v,_du),t.trian,t.quad)
end

function get_cell_residual(t::NonLinearFETerm,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  _uh = restrict(uh,t.trian)
  integrate(t.res(_uh,_v),t.trian,t.quad)
end

function get_cell_id(t::NonLinearFETerm)
  get_cell_id(t.trian)
end

"""
"""
struct AffineFETermFromCellMatVec <: AffineFETerm
  matvecfun::Function
  trian::Triangulation
end

function get_cell_matrix(t::AffineFETermFromCellMatVec,v,u)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  _v = restrict(v,t.trian)
  _u = restrict(u,t.trian)
  cellmatvec = t.matvecfun(_v,_u)
  cellmat, _ = unpair_arrays(cellmatvec)
  cellmat
end

function get_cell_vector(t::AffineFETermFromCellMatVec,v,uhd)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  trial = TrialFESpace(get_fe_space(uhd))
  u = get_cell_basis(trial)
  _v = restrict(v,t.trian)
  _u = restrict(u,t.trian)
  _cellvals = get_cell_values(t,uhd)
  cellmatvec = t.matvecfun(_v,_u)
  cellmatvec_with_diri = attach_dirichlet_bcs(cellmatvec,_cellvals)
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

function get_cell_values(t::AffineFETermFromCellMatVec,uhd)
  @assert is_a_fe_function(uhd)
  cellvals = get_cell_values(uhd)
  reindex(cellvals,t.trian)
end

function get_cell_matrix_and_vector(t::AffineFETermFromCellMatVec,v,u,uhd)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  _v = restrict(v,t.trian)
  _u = restrict(u,t.trian)
  _cellvals = get_cell_values(t,uhd)
  cellmatvec = t.matvecfun(_v,_u)
  cellmatvec_with_diri = attach_dirichlet_bcs(cellmatvec,_cellvals)
  (cellmatvec_with_diri, nothing, nothing)
end

struct FETermFromCellJacRes <: FETerm
  jacresfun::Function
  trian::Triangulation
end

function get_cell_jacobian(t::FETermFromCellJacRes,uh,v,du)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  _v = restrict(v,t.trian)
  _du = restrict(du,t.trian)
  _uh = restrict(uh,t.trian)
  celljacres = t.jacresfun(_uh,_v,_du)
  celljac, _ = unpair_arrays(celljacres)
  celljac
end

function get_cell_residual(t::FETermFromCellJacRes,uh,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  trial = TrialFESpace(get_fe_space(uh))
  du = get_cell_basis(trial)
  _v = restrict(v,t.trian)
  _du = restrict(du,t.trian)
  _uh = restrict(uh,t.trian)
  celljacres = t.jacresfun(_uh,_v,_du)
  _, cellres = unpair_arrays(celljacres)
  cellres
end

function get_cell_id(t::FETermFromCellJacRes)
  get_cell_id(t.trian)
end

function get_cell_jacobian_and_residual(t::FETermFromCellJacRes,uh,v,du)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  _v = restrict(v,t.trian)
  _du = restrict(du,t.trian)
  _uh = restrict(uh,t.trian)
  celljacres = t.jacresfun(_uh,_v,_du)
  (celljacres, nothing, nothing)
end

