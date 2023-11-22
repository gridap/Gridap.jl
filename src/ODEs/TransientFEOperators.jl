# TODO Conceptually, TransientFEOperator is a subtype of FEOperator and
# ODEOperator. Should there be a hierarchy of types?
"""
    abstract type TransientFEOperator <: GridapType end

Transient version of `FEOperator`
"""
abstract type TransientFEOperator{C<:ODEOperatorType} <: GridapType end

"""
    ODEOperatorType(::Type{<:TransientFEOperator})

Return the `ODEOperatorType` of the TransientFEOperator
"""
ODEOperatorType(::Type{<:TransientFEOperator{C}}) where {C} = C

"""
Return the assembler of the TransientFEOperator, which is constant for all time
steps.

Note: adaptive FE spaces involve generating new FE spaces and corresponding
operators, due to the immutable approach in `Gridap`
"""
get_assembler(op::TransientFEOperator) = @abstractmethod

# FEOperator interface
function FESpaces.get_test(op::TransientFEOperator)
  @abstractmethod
end

function FESpaces.get_trial(op::TransientFEOperator)
  @abstractmethod
end

"""
    get_algebraic_operator(op::TransientFEOperator)

Return an `ODEOperator` wrapper for the `TransientFEOperator` that can be
straightforwardly used with an `ODESolver`.
"""
function FESpaces.get_algebraic_operator(op::TransientFEOperator{C}) where {C}
  ODEOpFromFEOp{C}(op)
end

# ODEOperator interface
function Polynomials.get_order(op::TransientFEOperator)
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

function Algebra.allocate_residual(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  @abstractmethod
end

function Algebra.residual(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  r = allocate_residual(op, t, uh, cache)
  residual!(r, op, t, uh, cache)
  r
end

function Algebra.residual!(
  r::AbstractVector, op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  @abstractmethod
end

function Algebra.allocate_jacobian(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  cache
)
  @notimplemented
end

function Algebra.jacobian(
  op::TransientFEOperator,
  t::Real, uh::CellField,
  i::Int, γ::Real,
  cache
)
  J = allocate_jacobian(op, t, uh, cache)
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, op, t, uh, i, γ, cache)
  J
end

function Algebra.jacobian!(
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

###################################
# TransientFEOperatorFromWeakForm #
###################################
"""
    struct TransientFEOperatorFromWeakForm <: TransientFEOperator

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

# Default constructors
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

FESpaces.get_test(op::TransientFEOperatorFromWeakForm) = op.test

FESpaces.get_trial(op::TransientFEOperatorFromWeakForm) = op.trial

Polynomials.get_order(op::TransientFEOperatorFromWeakForm) = op.order

function Algebra.allocate_residual(
  op::TransientFEOperatorFromWeakForm,
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

function Algebra.residual!(
  r::AbstractVector, op::TransientFEOperatorFromWeakForm,
  t::Real, uh::CellField,
  cache
)
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V, op.res(t, uh, v))
  assemble_vector!(r, get_assembler(op), vecdata)
  r
end

function Algebra.allocate_jacobian(
  op::TransientFEOperatorFromWeakForm,
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

function Algebra.jacobian!(
  J::AbstractMatrix, op::TransientFEOperatorFromWeakForm,
  t::Real, uh::CellField,
  i::Int, γ::Real,
  cache
)
  matdata = _matdata_jacobian(op, t, uh, i, γ)
  assemble_matrix_add!(J, get_assembler(op), matdata)
  J
end

function jacobians!(
  J::AbstractMatrix, op::TransientFEOperatorFromWeakForm,
  t::Real, uh::CellField,
  γs::Tuple{Vararg{Real}},
  cache
)
  _matdata = ()
  for i in 0:get_order(op)
    γ = γs[i+1]
    if !iszero(γ)
      _matdata = (_matdata..., _matdata_jacobian(op, t, uh, i, γ))
    end
  end

  matdata = _vcat_matdata(_matdata)
  assemble_matrix_add!(J, get_assembler(op), matdata)
  J
end

#################
# ODEOpFromFEOp #
#################
"""
    struct ODEOpFromFEOp <: ODEOperator

Wrapper of `TransientFEOperator` that transforms it into a `ODEOperator`, i.e.
takes A(t, uh, ∂tuh, ..., ∂t^Nuh, vh) and returns A(t, uF, ∂tuF, ..., ∂t^NuF)
where uF, ∂tuF, ..., ∂t^NuF represent the free values of the
`EvaluationFunction`s uh, ∂tuh, ..., ∂t^Nuh
"""
struct ODEOpFromFEOp{C} <: ODEOperator{C}
  feop::TransientFEOperator{C}
end

# ODEOperator interface
Polynomials.get_order(op::ODEOpFromFEOp) = get_order(op.feop)

function allocate_cache(op::ODEOpFromFEOp)
  Ut = get_trial(op.feop)
  U = allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for i in 1:get_order(op)
    Uts = (Uts..., ∂t(Uts[i]))
    Us = (Us..., allocate_space(Uts[i+1]))
  end
  fe_cache = allocate_cache(op.feop)
  ode_cache = (Us, Uts, fe_cache)
  ode_cache
end

function allocate_cache(op::ODEOpFromFEOp, v::AbstractVector)
  cache = allocate_cache(op)
  _v = similar(v)
  (_v, cache)
end

function allocate_cache(op::ODEOpFromFEOp, v::AbstractVector, a::AbstractVector)
  cache = allocate_cache(op)
  _v = similar(v)
  _a = similar(a)
  (_v, _a, cache)
end

function update_cache!(cache, op::ODEOpFromFEOp, t::Real)
  _Us, Uts, fe_cache = cache
  Us = ()
  for i in 1:get_order(op)+1
    Us = (Us..., evaluate!(_Us[i], Uts[i], t))
  end
  fe_cache = update_cache!(fe_cache, op.feop, t)
  (Us, Uts, fe_cache)
end

function Algebra.allocate_residual(
  op::ODEOpFromFEOp,
  t::Real, u::AbstractVector,
  cache
)
  Us, Uts, fe_cache = cache
  uh = EvaluationFunction(Us[1], u)
  allocate_residual(op.feop, t, uh, fe_cache)
end

function Algebra.residual!(
  r::AbstractVector, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  residual!(r, op.feop, t, uh, cache)
end

function Algebra.allocate_jacobian(
  op::ODEOpFromFEOp,
  t::Real, u::AbstractVector,
  cache
)
  Us, Uts, fecache = cache
  uh = EvaluationFunction(Us[1], u)
  allocate_jacobian(op.feop, t, uh, fecache)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  i::Integer, γ::Real,
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  jacobian!(J, op.feop, t, uh, i, γ, cache)
end

function jacobians!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  jacobians!(J, op.feop, t, uh, γs, cache)
end

#########
# Utils #
#########
function _matdata_jacobian(
  op::TransientFEOperatorFromWeakForm,
  t::Real, uh::CellField,
  i::Integer, γ::Real
)
  Ut = evaluate(get_trial(op), nothing)
  V = get_test(op)
  du = get_trial_fe_basis(Ut)
  v = get_fe_basis(V)
  collect_cell_matrix(Ut, V, γ * op.jacs[i+1](t, uh, du, v))
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

  (term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
end

function _make_uh_from_us(op, us, cache)
  Ut, = cache

  u = EvaluationFunction(Ut[1], us[1])

  dus = ()
  for i in 2:get_order(op)+1
    dus = (dus..., EvaluationFunction(Ut[i], us[i]))
  end

  TransientCellField(u, dus)
end

########
# Test #
########
"""
    test_transient_fe_operator(op::TransientFEOperator, uh) -> Bool

Test the interface of `TransientFEOperator` specializations
"""
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

  jacobian!(J, op, t, xh, 0, 1, cache)
  @test J isa AbstractMatrix

  jacobian!(J, op, t, xh, 1, 1, cache)
  @test J isa AbstractMatrix

  jacobians!(J, op, t, xh, (1, 1), cache)
  @test J isa AbstractMatrix

  cache = update_cache!(cache, op, t)

  true
end
