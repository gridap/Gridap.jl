# QUESTION conceptually, TransientFEOperator is a subtype of FEOperator and
# ODEOperator. Should there be a hierarchy of types?
"""
A transient version of the `Gridap` `FEOperator` that depends on time
"""
abstract type TransientFEOperator{C<:ODEOperatorType} <: GridapType end

ODEOperatorType(::Type{<:TransientFEOperator{C}}) where {C} = C

"""
Return the assembler of the TransientFEOperator, which is constant for all time
steps.

Note: adaptive FE spaces involve generating new FE spaces and corresponding
operators, due to the immutable approach in `Gridap`
"""
get_assembler(op::TransientFEOperator) = @abstractmethod

# FEOperator interface
function get_test(op::TransientFEOperator)
  @abstractmethod
end

function get_trial(op::TransientFEOperator)
  @abstractmethod
end

"""
Return an `ODEOperator` wrapper for the `TransientFEOperator` that can be
straightforwardly used within `ODETools`.
"""
function get_algebraic_operator(op::TransientFEOperator{C}) where {C}
  ODEOpFromFEOp{C}(op)
end

# ODEOperator interface
function get_order(op::TransientFEOperator)
  @abstractmethod
end


function allocate_residual(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  @abstractmethod
end

function residual!(
  r::AbstractVector, op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  @abstractmethod
end

function allocate_jacobian(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  @notimplemented
end

function jacobian!(
  J::AbstractMatrix, op::TransientFEOperator,
  t::Real, uh::CellField,
  i::Int, γ::Real,
  cache
)
  @abstractmethod
end

function jacobians!(
  J::AbstractMatrix, op::TransientFEOperator,
  t::Real, uh::CellField,
  γs::Tuple{Vararg{Real}},
  cache
)
  @abstractmethod
end

# @fverdugo This function is just in case we need to override it in the future
# for some specialization. This default implementation is enough for the moment.
function allocate_cache(op::TransientFEOperator)
  nothing
end

function update_cache!(cache, op::TransientFEOperator, t::Real)
  nothing
end

###################################
# TransientFEOperatorFromWeakForm #
###################################
"""
Transient FE operator that is defined by a transient weak form
"""
struct TransientFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  res::Function
  jacs::Tuple{Vararg{Function}}
  assembler::Assembler
  trial
  test::FESpace
  order::Integer
end

function TransientFEOperator(res::Function, trial, test; order::Integer=1)
  function jac_0(t, u, du, dv)
    function res_0(y)
      u0 = TransientCellField(y, u.derivatives)
      res(t, u0, dv)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for i in 1:order
    function jac_i(t, u, dui, dv)
      function res_i(y)
        derivatives = (u.derivatives[1:i-1]..., y, u.derivatives[i+1:end]...)
        ui = TransientCellField(u.cellfield, derivatives)
        res(t, ui, dv)
      end
      jacobian(res_i, u.derivatives[i])
    end
    jacs = (jacs..., jac_i)
  end

  TransientFEOperator(res, jacs..., trial, test)
end

function TransientFEOperator(
  res::Function, jac::Function, jac_t::Function,
  trial, test
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientFEOperatorFromWeakForm{NonlinearODE}(
    res, (jac, jac_t),
    assembler, trial, test, 1
  )
end

function TransientFEOperator(
  res::Function, jac::Function, jac_t::Function, jac_tt::Function,
  trial, test
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientFEOperatorFromWeakForm{NonlinearODE}(
    res, (jac, jac_t, jac_tt),
    assembler, trial, test, 2
  )
end

# TransientFEOperator interface
get_assembler(op::TransientFEOperatorFromWeakForm) = op.assembler

get_test(op::TransientFEOperatorFromWeakForm) = op.test

get_trial(op::TransientFEOperatorFromWeakForm) = op.trial

get_order(op::TransientFEOperatorFromWeakForm) = op.order

function allocate_residual(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  duhs = ()
  for i in 1:get_order(op)
    duhs = (duhs..., uh)
  end
  us = TransientCellField(uh, duhs)

  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V, op.res(t, us, v))
  allocate_vector(get_assembler(op), vecdata)
end

function residual!(
  r::AbstractVector, op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V, op.res(t, uh, v))
  assemble_vector!(r, get_assembler(op), vecdata)
  r
end

function allocate_jacobian(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  duhs = ()
  for i in 1:get_order(op)
    duhs = (duhs..., uh)
  end
  us = TransientCellField(uh, duhs)

  _matdata = ()
  for i in 0:get_order(op)
    _matdata = (_matdata..., _matdata_jacobian(op, t, us, i, 0))
  end

  matdata = _vcat_matdata(_matdata)
  allocate_matrix(get_assembler(op), matdata)
end

function jacobian!(
  J::AbstractMatrix, op::TransientFEOperator,
  t::Real, uh::CellField,
  i::Int, γ::Real,
  cache
)
  matdata = _matdata_jacobian(op, t, uh, i, γ)
  assemble_matrix_add!(J, get_assembler(op), matdata)
  J
end

function jacobians!(
  J::AbstractMatrix, op::TransientFEOperator,
  t::Real, uh::CellField,
  γs::Tuple{Vararg{Real}},
  cache
)
  _matdata = ()
  for i in 0:get_order(op)
    if !iszero(γ[i])
      _matdata = (_matdata..., _matdata_jacobian(op, t, uh, i, γs[i]))
    end
  end

  matdata = _vcat_matdata(_matdata)
  assemble_matrix_add!(J, get_assembler(op), matdata)
  J
end

#########
# Utils #
#########
function _matdata_jacobian(
  op::TransientFEOperatorsFromWeakForm,
  t::Real, uh::CellField,
  i::Integer, γ::Real
)
  Ut = evaluate(get_trial(op), nothing)
  V = get_test(op)
  du = get_trial_fe_basis(Ut)
  v = get_fe_basis(V)
  matdata = collect_cell_matrix(Ut, V, γ * op.jacs[i+1](t, uh, du, v))
end

function _vcat_matdata(_matdata)
  term_to_cellmat_j = ()
  term_to_cellidsrows_j = ()
  term_to_cellidscols_j = ()
  for j in 1:length(_matdata)
    term_to_cellmat_j = (term_to_cellmat_j..., _matdata[j][1])
    term_to_cellidsrows_j = (term_to_cellidsrows_j..., _matdata[j][2])
    term_to_cellidscols_j = (term_to_cellidscols_j..., _matdata[j][3])
  end

  term_to_cellmat = vcat(term_to_cellmat_j...)
  term_to_cellidsrows = vcat(term_to_cellidsrows_j...)
  term_to_cellidscols = vcat(term_to_cellidscols_j...)

  matdata = (term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
end

########
# Test #
########
function test_transient_fe_operator(op::TransientFEOperator, uh)
  odeop = get_algebraic_operator(op)
  @test odeop isa ODEOperator

  V = get_test(op)
  @test V isa FESpace

  t = 0.0
  U = get_trial(op)
  Ut = U(t)
  @test isa(Ut, FESpace)

  cache = allocate_cache(op)

  r = allocate_residual(op, t, uh, cache)
  @test r isa AbstractVector

  xh = TransientCellField(uh, (uh,))

  residual!(r, op, t, xh, cache)
  @test r isa AbstractVector

  J = allocate_jacobian(op, t, uh, cache)
  @test J isa AbstractMatrix

  jacobian!(J, op, t, xh, 1, 1, cache)
  @test J isa AbstractMatrix

  jacobian!(J, op, t, xh, 2, 1, cache)
  @test J isa AbstractMatrix

  jacobians!(J, op, t, xh, (1, 1), cache)
  @test J isa AbstractMatrix

  cache = update_cache!(cache, op, t)

  true
end
