"""
A transient version of the `Gridap` `FEOperator` that depends on time
"""
abstract type TransientFEOperator{C<:OperatorType} <: GridapType end

"""
Returns the test space
"""
function get_test(op::TransientFEOperator)
  @abstractmethod
end

"""
Returns the (possibly) time-dependent trial space
"""
function get_trial(op::TransientFEOperator)
  @abstractmethod # time dependent
end

function allocate_residual(op::TransientFEOperator,t0,uh,cache)
  @abstractmethod
end

function allocate_jacobian(op::TransientFEOperator,t0,uh,cache)
  @notimplemented
end

"""
Idem as `residual!` of `ODEOperator`
"""
function residual!(
  b::AbstractVector,
  op::TransientFEOperator,
  t::Real,
  xh::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  cache)
  @abstractmethod
end

"""
Idem as `jacobian!` of `ODEOperator`
"""
function jacobian!(
  A::AbstractMatrix,
  op::TransientFEOperator,
  t::Real,
  xh::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  i::Int,
  γᵢ::Real,
  cache)
  @abstractmethod
end

"""
Idem as `jacobians!` of `ODEOperator`
"""
function jacobians!(
  A::AbstractMatrix,
  op::TransientFEOperator,
  t::Real,
  x::Tuple{Vararg{AbstractVector}},
  γ::Tuple{Vararg{Real}},
  cache)
  @abstractmethod
end

"""
Returns the assembler, which is constant for all time steps for a given FE
operator.

Note: adaptive FE spaces involve to generate new FE spaces and
corresponding operators, due to the ummutable approach in `Gridap`
"""
get_assembler(feop::TransientFEOperator) = @abstractmethod


# Default API

"""
Returns a `ODEOperator` wrapper of the `TransientFEOperator` that can be
straightforwardly used with the `ODETools` module.
"""
function get_algebraic_operator(feop::TransientFEOperator{C}) where C
  ODEOpFromFEOp{C}(feop)
end

OperatorType(::Type{<:TransientFEOperator{C}}) where C = C

# @fverdugo This function is just in case we need to override it in the future for some specialization.
# This default implementation is enough for the moment.
function allocate_cache(op::TransientFEOperator)
  nothing
end

function update_cache!(cache::Nothing,op::TransientFEOperator,t::Real)
  nothing
end

# Specializations

"""
Transient FE operator that is defined by a transient Weak form
"""
struct TransientFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  res::Function
  rhs::Function
  jacs::Tuple{Vararg{Function}}
  assem_t::Assembler
  trials::Tuple{Vararg{Any}}
  test::FESpace
  order::Integer
end

function TransientConstantFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,u,v) = -1.0 * b(v)
  rhs(t,u,v) = b(v)
  jac(t,u,du,v) = a(du,v)
  jac_t(t,u,dut,v) = m(dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Constant}(res,rhs,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientConstantMatrixFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,u,v) = m(∂t(u),v) + a(u,v) - b(t,v)
  rhs(t,u,v) = b(t,v) - a(u,v)
  jac(t,u,du,v) = a(du,v)
  jac_t(t,u,dut,v) = m(dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{ConstantMatrix}(res,rhs,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientAffineFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,u,v) = m(t,∂t(u),v) + a(t,u,v) - b(t,v)
  rhs(t,u,v) = b(t,v) - a(t,u,v)
  jac(t,u,du,v) = a(t,du,v)
  jac_t(t,u,dut,v) = m(t,dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Affine}(res,rhs,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Nonlinear}(res,rhs_error,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end


function TransientConstantFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,u,v) = -1.0 * b(v)
  rhs(t,u,v) = b(v)
  jac(t,u,du,v) = a(du,v)
  jac_t(t,u,dut,v) = c(dut,v)
  jac_tt(t,u,dutt,v) = m(dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{Constant}(
    res,rhs,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientConstantMatrixFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,u,v) = m(∂tt(u),v) + c(∂t(u),v) + a(u,v) - b(t,v)
  rhs(t,u,v) = b(t,v) - c(∂t(u),v) - a(u,v)
  jac(t,u,du,v) = a(du,v)
  jac_t(t,u,dut,v) = c(dut,v)
  jac_tt(t,u,dutt,v) = m(dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{ConstantMatrix}(
    res,rhs,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientAffineFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,u,v) = m(t,∂tt(u),v) + c(t,∂t(u),v) + a(t,u,v) - b(t,v)
  rhs(t,u,v) = b(t,v) - c(t,∂t(u),v) - a(t,u,v)
  jac(t,u,du,v) = a(t,du,v)
  jac_t(t,u,dut,v) = c(t,dut,v)
  jac_tt(t,u,dutt,v) = m(t,dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{Affine}(
    res,rhs,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  jac_tt::Function,trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{Nonlinear}(
    res,rhs_error,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientFEOperator(res::Function,trial,test;order::Integer=1)
  function jac_0(t,x,dx0,dv)
    function res_0(y)
      x0 = TransientCellField(y,x.derivatives)
      res(t,x0,dv)
    end
    jacobian(res_0,x.cellfield)
  end
  jacs = (jac_0,)
  for i in 1:order
    function jac_i(t,x,dxi,dv)
      function res_i(y)
        derivatives = (x.derivatives[1:i-1]...,y,x.derivatives[i+1:end]...)
        xi = TransientCellField(x.cellfield,derivatives)
        res(t,xi,dv)
      end
      jacobian(res_i,x.derivatives[i])
    end
    jacs = (jacs...,jac_i)
  end
  TransientFEOperator(res,jacs...,trial,test)
end

function allocate_residual(
  op::TransientFEOperatorFromWeakForm,
  t0::Real,
  uh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  dxh = ()
  for i in 1:get_order(op)
    dxh = (dxh...,uh)
  end
  xh = TransientCellField(uh,dxh)
  vecdata = collect_cell_vector(V,op.res(t0,xh,v))
  allocate_vector(op.assem_t,vecdata)
end

function residual!(
  b::AbstractVector,
  op::TransientFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.res(t,xh,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

"""
Transient FE operator that is defined by a transient Weak form with the
form: LHS(t,u,∂u/∂t,...) ∂u/∂t = RHS(t,u,∂u/∂t,...). Used in Runge-Kutta schemes
"""
struct TransientRKFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  lhs::Function
  rhs::Function
  jacs::Tuple{Vararg{Function}}
  assem_t::Assembler
  trials::Tuple{Vararg{Any}}
  test::FESpace
  order::Integer
end

function TransientRungeKuttaFEOperator(lhs::Function,rhs::Function,jac::Function,
  jac_t::Function,trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientRKFEOperatorFromWeakForm{Nonlinear}(lhs,rhs,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientRungeKuttaFEOperator(lhs::Function,rhs::Function,trial,test)
  res(t,u,v) = lhs(t,u,v) - rhs(t,u,v)
  function jac_0(t,x,dx0,dv)
    function res_0(y)
      x0 = TransientCellField(y,x.derivatives)
      res(t,x0,dv)
    end
    jacobian(res_0,x.cellfield)
  end
  jacs = (jac_0,)
  function jac_t(t,x,dxt,dv)
    function res_t(y)
      derivatives = (y,x.derivatives[2:end]...)
      xt = TransientCellField(x.cellfield,derivatives)
      res(t,xt,dv)
    end
    jacobian(res_t,x.derivatives[1])
  end
  jacs = (jac_0,jac_t)
  TransientRungeKuttaFEOperator(lhs,rhs,jacs...,trial,test)
end

function allocate_residual(
  op::TransientRKFEOperatorFromWeakForm,
  t0::Real,
  uh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  dxh = ()
  for i in 1:get_order(op)
    dxh = (dxh...,uh)
  end
  xh = TransientCellField(uh,dxh)
  vecdata = collect_cell_vector(V,op.lhs(t0,xh,v))
  allocate_vector(op.assem_t,vecdata)
end

function lhs!(
  b::AbstractVector,
  op::TransientRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.lhs(t,xh,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function rhs!(
  rhs::AbstractVector,
  op::TransientRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.rhs(t,xh,v))
  assemble_vector!(rhs,op.assem_t,vecdata)
  rhs
end

# IMEX-RK Transient FE operators
"""
Transient FE operator that is defined by a transient Weak form with the
form: LHS(t,u,∂u/∂t,...) ∂u/∂t = I_RHS(t,u,∂u/∂t,...) + E_RHS(t,u,∂u/∂t,...).
Used in Implicit-Explicit Runge-Kutta schemes
"""
struct TransientIMEXRKFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  lhs::Function
  rhs::Function
  explicit_rhs::Function
  jacs::Tuple{Vararg{Function}}
  assem_t::Assembler
  trials::Tuple{Vararg{Any}}
  test::FESpace
  order::Integer
end

function TransientIMEXRungeKuttaFEOperator(lhs::Function,rhs::Function,
  explicit_rhs::Function,jac::Function,jac_t::Function,trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientIMEXRKFEOperatorFromWeakForm{Nonlinear}(lhs,rhs,explicit_rhs,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientIMEXRungeKuttaFEOperator(lhs::Function,rhs::Function,
  explicit_rhs::Function,trial,test)
  res(t,u,v) = lhs(t,u,v) - rhs(t,u,v)
  function jac_0(t,x,dx0,dv)
    function res_0(y)
      x0 = TransientCellField(y,x.derivatives)
      res(t,x0,dv)
    end
    jacobian(res_0,x.cellfield)
  end
  jacs = (jac_0,)
  function jac_t(t,x,dxt,dv)
    function res_t(y)
      derivatives = (y,x.derivatives[2:end]...)
      xt = TransientCellField(x.cellfield,derivatives)
      res(t,xt,dv)
    end
    jacobian(res_t,x.derivatives[1])
  end
  jacs = (jac_0,jac_t)
  TransientIMEXRungeKuttaFEOperator(lhs,rhs,explicit_rhs,jacs...,trial,test)
end

function allocate_residual(
  op::TransientIMEXRKFEOperatorFromWeakForm,
  t0::Real,
  uh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  dxh = ()
  for i in 1:get_order(op)
    dxh = (dxh...,uh)
  end
  xh = TransientCellField(uh,dxh)
  vecdata = collect_cell_vector(V,op.lhs(t0,xh,v))
  allocate_vector(op.assem_t,vecdata)
end

function lhs!(
  b::AbstractVector,
  op::TransientIMEXRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.lhs(t,xh,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function rhs!(
  rhs::AbstractVector,
  op::TransientIMEXRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.rhs(t,xh,v))
  assemble_vector!(rhs,op.assem_t,vecdata)
  rhs
end

function explicit_rhs!(
  explicit_rhs::AbstractVector,
  op::TransientIMEXRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.explicit_rhs(t,xh,v))
  assemble_vector!(explicit_rhs,op.assem_t,vecdata)
  explicit_rhs
end


# EX-RK Transient FE operators
"""
Transient FE operator that is defined by a transient Weak form with the
form: LHS(t,u,∂u/∂t,...) ∂u/∂t = E_RHS(t,u,∂u/∂t,...).
Used in Explicit Runge-Kutta schemes
"""
struct TransientEXRKFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  lhs::Function
  rhs::Function
  jacs::Tuple{Vararg{Function}}
  assem_t::Assembler
  trials::Tuple{Vararg{Any}}
  test::FESpace
  order::Integer
end

function TransientEXRungeKuttaFEOperator(lhs::Function,
  rhs::Function,jac::Function,jac_t::Function,trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientEXRKFEOperatorFromWeakForm{Nonlinear}(lhs,rhs,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientEXRungeKuttaFEOperator(lhs::Function,
  rhs::Function,trial,test)
  res(t,u,v) = lhs(t,u,v) - rhs(t,u,v)
  function jac_0(t,x,dx0,dv)
    function res_0(y)
      x0 = TransientCellField(y,x.derivatives)
      res(t,x0,dv)
    end
    jacobian(res_0,x.cellfield)
  end
  jacs = (jac_0,)
  function jac_t(t,x,dxt,dv)
    function res_t(y)
      derivatives = (y,x.derivatives[2:end]...)
      xt = TransientCellField(x.cellfield,derivatives)
      res(t,xt,dv)
    end
    jacobian(res_t,x.derivatives[1])
  end
  jacs = (jac_0,jac_t)
  TransientEXRungeKuttaFEOperator(lhs,rhs,jacs...,trial,test)
end

function allocate_residual(
  op::TransientEXRKFEOperatorFromWeakForm,
  t0::Real,
  uh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  dxh = ()
  for i in 1:get_order(op)
    dxh = (dxh...,uh)
  end
  xh = TransientCellField(uh,dxh)
  vecdata = collect_cell_vector(V,op.lhs(t0,xh,v))
  allocate_vector(op.assem_t,vecdata)
end

function lhs!(
  b::AbstractVector,
  op::TransientEXRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.lhs(t,xh,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function rhs!(
  rhs::AbstractVector,
  op::TransientEXRKFEOperatorFromWeakForm,
  t::Real,
  xh::T,
  cache) where T
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V,op.rhs(t,xh,v))
  assemble_vector!(rhs,op.assem_t,vecdata)
  rhs
end

# function explicit_rhs!(
#   explicit_rhs::AbstractVector,
#   op::TransientEXRKFEOperatorFromWeakForm,
#   t::Real,
#   xh::T,
#   cache) where T
#   V = get_test(op)
#   v = get_fe_basis(V)
#   vecdata = collect_cell_vector(V,op.explicit_rhs(t,xh,v))
#   assemble_vector!(explicit_rhs,op.assem_t,vecdata)
#   explicit_rhs
# end








# Common functions

TransientFEOperatorsFromWeakForm = Union{TransientFEOperatorFromWeakForm,
TransientRKFEOperatorFromWeakForm, TransientIMEXRKFEOperatorFromWeakForm,
TransientEXRKFEOperatorFromWeakForm}

function SparseMatrixAssembler(
  trial::Union{TransientTrialFESpace,TransientMultiFieldTrialFESpace},
  test::FESpace)
  SparseMatrixAssembler(evaluate(trial,nothing),test)
end

get_assembler(op::TransientFEOperatorsFromWeakForm) = op.assem_t

get_test(op::TransientFEOperatorsFromWeakForm) = op.test

get_trial(op::TransientFEOperatorsFromWeakForm) = op.trials[1]

get_order(op::TransientFEOperatorsFromWeakForm) = op.order


function allocate_jacobian(
  op::TransientFEOperatorsFromWeakForm,
  t0::Real,
  uh::CellField,
  cache)
  _matdata_jacobians = fill_initial_jacobians(op,t0,uh)
  matdata = _vcat_matdata(_matdata_jacobians)
  allocate_matrix(op.assem_t,matdata)
end

function jacobian!(
  A::AbstractMatrix,
  op::TransientFEOperatorsFromWeakForm,
  t::Real,
  xh::T,
  i::Integer,
  γᵢ::Real,
  cache) where T
  matdata = _matdata_jacobian(op,t,xh,i,γᵢ)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobians!(
  A::AbstractMatrix,
  op::TransientFEOperatorsFromWeakForm,
  t::Real,
  xh::TransientCellField,
  γ::Tuple{Vararg{Real}},
  cache)
  _matdata_jacobians = fill_jacobians(op,t,xh,γ)
  matdata = _vcat_matdata(_matdata_jacobians)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function fill_initial_jacobians(op::TransientFEOperatorsFromWeakForm,t0::Real,uh)
  dxh = ()
  for i in 1:get_order(op)
    dxh = (dxh...,uh)
  end
  xh = TransientCellField(uh,dxh)
  _matdata = ()
  for i in 1:get_order(op)+1
    _matdata = (_matdata...,_matdata_jacobian(op,t0,xh,i,0.0))
  end
  return _matdata
end

function fill_jacobians(
  op::TransientFEOperatorsFromWeakForm,
  t::Real,
  xh::T,
  γ::Tuple{Vararg{Real}}) where T
  _matdata = ()
  for i in 1:get_order(op)+1
    if (γ[i] > 0.0)
      _matdata = (_matdata...,_matdata_jacobian(op,t,xh,i,γ[i]))
    end
  end
  return _matdata
end

function _vcat_matdata(_matdata)
  term_to_cellmat_j = ()
  term_to_cellidsrows_j = ()
  term_to_cellidscols_j = ()
  for j in 1:length(_matdata)
    term_to_cellmat_j = (term_to_cellmat_j...,_matdata[j][1])
    term_to_cellidsrows_j = (term_to_cellidsrows_j...,_matdata[j][2])
    term_to_cellidscols_j = (term_to_cellidscols_j...,_matdata[j][3])
  end

  term_to_cellmat = vcat(term_to_cellmat_j...)
  term_to_cellidsrows = vcat(term_to_cellidsrows_j...)
  term_to_cellidscols = vcat(term_to_cellidscols_j...)

  matdata = (term_to_cellmat,term_to_cellidsrows, term_to_cellidscols)
end

function _matdata_jacobian(
  op::TransientFEOperatorsFromWeakForm,
  t::Real,
  xh::T,
  i::Integer,
  γᵢ::Real) where T
  Uh = evaluate(get_trial(op),nothing)
  V = get_test(op)
  du = get_trial_fe_basis(Uh)
  v = get_fe_basis(V)
  matdata = collect_cell_matrix(Uh,V,γᵢ*op.jacs[i](t,xh,du,v))
end

function rhs_error(t::Real,xh,v)
  error("The \"rhs\" function is not defined for this TransientFEOperator.
  Please, try to use another type of TransientFEOperator that supports this
  functionality.")
end

# Tester

function test_transient_fe_operator(op::TransientFEOperator,uh)
  odeop = get_algebraic_operator(op)
  @test isa(odeop,ODEOperator)
  cache = allocate_cache(op)
  V = get_test(op)
  @test isa(V,FESpace)
  U = get_trial(op)
  U0 = U(0.0)
  @test isa(U0,FESpace)
  r = allocate_residual(op,0.0,uh,cache)
  @test isa(r,AbstractVector)
  xh = TransientCellField(uh,(uh,))
  residual!(r,op,0.0,xh,cache)
  @test isa(r,AbstractVector)
  J = allocate_jacobian(op,0.0,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,xh,1,1.0,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,xh,2,1.0,cache)
  @test isa(J,AbstractMatrix)
  jacobians!(J,op,0.0,xh,(1.0,1.0),cache)
  @test isa(J,AbstractMatrix)
  cache = update_cache!(cache,op,0.0)
  true
end
