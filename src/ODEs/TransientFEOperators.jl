# TODO Conceptually, TransientFEOperator is a subtype of FEOperator and
# ODEOperator. Should there be a hierarchy of types?
"""
    abstract type TransientFEOperator <: GridapType end

Transient version of `FEOperator` corresponding to a residual of the form
`residual(t, u, v) = 0`, that can involve time derivatives of `u`.
"""
abstract type TransientFEOperator{C<:ODEOperatorType} <: GridapType end

"""
    ODEOperatorType(::Type{<:TransientFEOperator}) -> ODEOperatorType

Return the `ODEOperatorType` of the `TransientFEOperator`
"""
ODEOperatorType(::Type{<:TransientFEOperator{C}}) where {C} = C

"""
    get_assembler(op::TransientFEOperator) -> Assembler

Return the assembler of the `TransientFEOperator`
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
    get_algebraic_operator(op::TransientFEOperator) -> ODEOperator

Return an `ODEOperator` wrapper around the `TransientFEOperator`, that can be
straightforwardly used with an `ODESolver`
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
  cache
end

# TODO replace CellField by TransientCellField?
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

#############################
# TransientFEOpFromWeakForm #
#############################
"""
    struct TransientFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
"""
struct TransientFEOpFromWeakForm{C} <: TransientFEOperator{C}
  res::Function
  jacs::Tuple{Vararg{Function}}
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Default constructors
function TransientFEOperator(res::Function, trial, test; order::Integer=1)
  if trial isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end

  function jac_0(t, u, du, dv)
    function res_0(y)
      u0 = TransientCellFieldType(y, u.derivatives)
      res(t, u0, dv)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for i in 1:order
    function jac_i(t, u, dui, dv)
      function res_i(y)
        derivatives = (u.derivatives[1:i-1]..., y, u.derivatives[i+1:end]...)
        ui = TransientCellFieldType(u.cellfield, derivatives)
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
  TransientFEOpFromWeakForm{NonlinearODE}(
    res, (jac, jac_t),
    assembler, trial, test, 1
  )
end

function TransientFEOperator(
  res::Function, jac::Function, jac_t::Function, jac_tt::Function,
  trial, test
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientFEOpFromWeakForm{NonlinearODE}(
    res, (jac, jac_t, jac_tt),
    assembler, trial, test, 2
  )
end

# TransientFEOperator interface
get_assembler(op::TransientFEOpFromWeakForm) = op.assembler

FESpaces.get_test(op::TransientFEOpFromWeakForm) = op.test

FESpaces.get_trial(op::TransientFEOpFromWeakForm) = op.trial

Polynomials.get_order(op::TransientFEOpFromWeakForm) = op.order

function Algebra.allocate_residual(
  op::TransientFEOpFromWeakForm,
  t::Real, uh::CellField,
  cache
)
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V, op.res(t, uh, v))
  allocate_vector(get_assembler(op), vecdata)
end

function Algebra.residual!(
  r::AbstractVector, op::TransientFEOpFromWeakForm,
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
  op::TransientFEOpFromWeakForm,
  t::Real, uh::CellField,
  cache
)
  _matdata = ()
  for i in 0:get_order(op)
    _matdata = (_matdata..., _matdata_jacobian(op, t, uh, i, 0))
  end

  matdata = _vcat_matdata(_matdata)
  allocate_matrix(get_assembler(op), matdata)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::TransientFEOpFromWeakForm,
  t::Real, uh::CellField,
  i::Int, γ::Real,
  cache
)
  matdata = _matdata_jacobian(op, t, uh, i, γ)
  assemble_matrix_add!(J, get_assembler(op), matdata)
  J
end

function jacobians!(
  J::AbstractMatrix, op::TransientFEOpFromWeakForm,
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

#######################################
# TransientMassLinearFEOpFromWeakForm #
#######################################
# TODO can also provide the mass bilinear form in terms of `u` and call it
# with the proper derivative order within `residual` and `jacobian`. The
# advantage is that the jacobian w.r.t. the highest-order time derivative
# does not need automatic differentiation
"""
    struct TransientMassLinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```math
mass(t, u, v) + res(t, u, v) = 0
```, where `mass` is linear in `∂t^N(u)` (`N` is the order of the operator) and
`res` does not depend on `∂t^N(u)` (`mass` and `res` can depend on lower-order
time-derivatives of `u` by applying the `∂t` operator)
"""
struct TransientMassLinearFEOpFromWeakForm{C<:AbstractMassLinearODE} <: TransientFEOperator{C}
  mass::Function
  res::Function
  jacs::Tuple{Vararg{Function}}
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Default constructors
function TransientMassLinearFEOperator(
  mass::Function, res::Function,
  trial, test; order::Integer=1
)
  if trial isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end

  function jac_0(t, u, du, dv)
    function res_0(y)
      u0 = TransientCellFieldType(y, u.derivatives)
      mass(t, u0, dv) + res(t, u0, dv)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for i in 1:order-1
    function jac_i(t, u, dui, dv)
      function res_i(y)
        derivatives = (u.derivatives[1:i-1]..., y, u.derivatives[i+1:end]...)
        ui = TransientCellFieldType(u.cellfield, derivatives)
        mass(t, ui, dv) + res(t, ui, dv)
      end
      jacobian(res_i, u.derivatives[i])
    end
    jacs = (jacs..., jac_i)
  end

  function jac_N(t, u, duN, dv)
    function res_N(y)
      derivatives = (u.derivatives[1:order-1]..., y)
      uN = TransientCellFieldType(u.cellfield, derivatives)
      mass(t, uN, dv)
    end
    jacobian(res_N, u.derivatives[order])
  end
  jacs = (jacs..., jac_N)

  TransientMassLinearFEOperator(mass, res, jacs..., trial, test)
end

function TransientMassLinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function,
  trial, test
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientMassLinearFEOpFromWeakForm{MassLinearODE}(
    mass, res, (jac, jac_t),
    assembler, trial, test, 1
  )
end

function TransientMassLinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientMassLinearFEOpFromWeakForm{MassLinearODE}(
    mass, res, (jac, jac_t, jac_tt),
    assembler, trial, test, 2
  )
end

# TransientFEOperator interface
get_assembler(op::TransientMassLinearFEOpFromWeakForm) = op.assembler

FESpaces.get_test(op::TransientMassLinearFEOpFromWeakForm) = op.test

FESpaces.get_trial(op::TransientMassLinearFEOpFromWeakForm) = op.trial

Polynomials.get_order(op::TransientMassLinearFEOpFromWeakForm) = op.order

function Algebra.allocate_residual(
  op::TransientMassLinearFEOpFromWeakForm,
  t::Real, uh::CellField,
  cache
)
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V, op.res(t, uh, v))
  # Could also take op.mass or op.mass + op.res above
  allocate_vector(get_assembler(op), vecdata)
end

function Algebra.residual!(
  r::AbstractVector, op::TransientMassLinearFEOpFromWeakForm,
  t::Real, uh::CellField,
  cache
)
  V = get_test(op)
  v = get_fe_basis(V)
  vecdata = collect_cell_vector(V, op.res(t, uh, v))
  assemble_vector!(r, get_assembler(op), vecdata)
  # Special case for the highest-order derivative
  vecdata = collect_cell_vector(V, op.mass(t, uh, v))
  assemble_vector_add!(r, get_assembler(op), vecdata)
  r
end

function Algebra.allocate_jacobian(
  op::TransientMassLinearFEOpFromWeakForm,
  t::Real, uh::CellField,
  cache
)
  _matdata = ()
  for i in 0:get_order(op)
    _matdata = (_matdata..., _matdata_jacobian(op, t, uh, i, 0))
  end

  matdata = _vcat_matdata(_matdata)
  allocate_matrix(get_assembler(op), matdata)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::TransientMassLinearFEOpFromWeakForm,
  t::Real, uh::CellField,
  i::Int, γ::Real,
  cache
)
  matdata = _matdata_jacobian(op, t, uh, i, γ)
  assemble_matrix_add!(J, get_assembler(op), matdata)
  J
end

function jacobians!(
  J::AbstractMatrix, op::TransientMassLinearFEOpFromWeakForm,
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
    struct ODEOpFromFEOp <: ODEOperator end

Wrapper that transforms a `TransientFEOperator` into an `ODEOperator`, i.e.
takes `residual(t, uh, ∂t(uh), ..., ∂t^N(uh), vh)` and returns
`residual(t, uf, ∂t(uf), ..., ∂t^N(uf))`, where `uf` represent the free values
of the `EvaluationFunction` uh
"""
struct ODEOpFromFEOp{C} <: ODEOperator{C}
  fe_op::TransientFEOperator{C}
end

# ODEOperator interface
Polynomials.get_order(op::ODEOpFromFEOp) = get_order(op.fe_op)

function allocate_cache(op::ODEOpFromFEOp)
  Ut = get_trial(op.fe_op)
  U = allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for i in 1:get_order(op)
    Uts = (Uts..., ∂t(Uts[i]))
    Us = (Us..., allocate_space(Uts[i+1]))
  end
  fe_cache = allocate_cache(op.fe_op)
  ode_cache = (Us, Uts, fe_cache)
  ode_cache
end

function update_cache!(cache, op::ODEOpFromFEOp, t::Real)
  _Us, Uts, fe_cache = cache
  Us = ()
  for i in 0:get_order(op)
    Us = (Us..., evaluate!(_Us[i+1], Uts[i+1], t))
  end
  fe_cache = update_cache!(fe_cache, op.fe_op, t)
  (Us, Uts, fe_cache)
end

function Algebra.allocate_residual(
  op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  _, _, fe_cache = cache
  uh = _make_uh_from_us(op.fe_op, us, cache)
  allocate_residual(op.fe_op, t, uh, fe_cache)
end

function Algebra.residual!(
  r::AbstractVector, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  residual!(r, op.fe_op, t, uh, cache)
end

function Algebra.allocate_jacobian(
  op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  _, _, fe_cache = cache
  uh = _make_uh_from_us(op.fe_op, us, cache)
  allocate_jacobian(op.fe_op, t, uh, fe_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  i::Integer, γ::Real,
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  jacobian!(J, op.fe_op, t, uh, i, γ, cache)
end

function jacobians!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  jacobians!(J, op.fe_op, t, uh, γs, cache)
end

#########
# Utils #
#########
function _matdata_jacobian(
  op::Union{TransientFEOpFromWeakForm,TransientMassLinearFEOpFromWeakForm},
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
  for i in 1:get_order(op)
    dus = (dus..., EvaluationFunction(Ut[i+1], us[i+1]))
  end

  if first(Ut) isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end
  TransientCellFieldType(u, dus)
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

  t = 0.0
  U = get_trial(op)
  Ut = U(t)
  V = get_test(op)

  @test V isa FESpace
  @test Ut isa FESpace

  if U isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end

  uhₜ = TransientCellFieldType(uh, (uh,))
  cache = allocate_cache(op)

  r = allocate_residual(op, t, uhₜ, cache)
  @test r isa AbstractVector

  residual!(r, op, t, uhₜ, cache)
  @test r isa AbstractVector

  J = allocate_jacobian(op, t, uhₜ, cache)
  @test J isa AbstractMatrix

  jacobian!(J, op, t, uhₜ, 0, 1, cache)
  @test J isa AbstractMatrix

  jacobian!(J, op, t, uhₜ, 1, 1, cache)
  @test J isa AbstractMatrix

  jacobians!(J, op, t, uhₜ, (1, 1), cache)
  @test J isa AbstractMatrix

  cache = update_cache!(cache, op, t)

  true
end
