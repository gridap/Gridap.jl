"""
    abstract type TransientFEOperator <: GridapType end

Transient version of `FEOperator` corresponding to a residual of the form
`residual(t, u, v) = 0`, that can involve time derivatives of `u`.

Important note: For now, the residual and jacobians cannot be directly computed
on a `TransientFEOperator`. They have to be evaluated on the corresponding
algebraic operator, which is an `ODEOperator`. As such, `TransientFEOperator` is
not exactly a subtype of `FEOperator`, but rather at the intersection of
`FEOperator` and `ODEOperator`.

# Mandatory
- [`get_test(op)`](@ref)
- [`get_trial(op)`](@ref)
- [`get_algebraic_operator(op)`](@ref)
- [`get_order(op)`](@ref)
- [`get_assembler(op)`](@ref)
- [`get_res(op::TransientFEOperator)`](@ref)
- [`get_jacs(op::TransientFEOperator)`](@ref)
- [`get_mass(op::TransientFEOperator{MassLinearODE})`](@ref)
- [`get_forms(op::TransientFEOperator{LinearODE})`](@ref)

# Optional
- [`is_jacobian_constant(op, k)`](@ref)
- [`is_forcing_constant(op::TransientFEOperator{<:AbstractMassLinearODE})`](@ref)
"""
abstract type TransientFEOperator{C<:ODEOperatorType} <: GridapType end

"""
    ODEOperatorType(::Type{<:TransientFEOperator}) -> ODEOperatorType

Return the `ODEOperatorType` of the `TransientFEOperator`.
"""
ODEOperatorType(op::TransientFEOperator) = ODEOperatorType(typeof(op))
ODEOperatorType(::Type{<:TransientFEOperator{C}}) where {C} = C

# FEOperator interface
function FESpaces.get_test(op::TransientFEOperator)
  @abstractmethod
end

function FESpaces.get_trial(op::TransientFEOperator)
  @abstractmethod
end

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

function is_jacobian_constant(op::TransientFEOperator, k::Integer)
  false
end

function is_forcing_constant(op::TransientFEOperator{<:AbstractMassLinearODE})
  false
end

# TransientFEOperator interface
"""
    get_assembler(op::TransientFEOperator) -> Assembler

Return the assembler of the `TransientFEOperator`.
"""
function get_assembler(op::TransientFEOperator)
  @abstractmethod
end

"""
    get_res(op::TransientFEOperator) -> Function

Return the (partial) residual of the `TransientFEOperator`.
"""
function get_res(op::TransientFEOperator)
  @abstractmethod
end

"""
    get_jacs(op::TransientFEOperator) -> Tuple{Vararg{Function}}

Return the jacobians of the `TransientFEOperator`.
"""
function get_jacs(op::TransientFEOperator)
  @abstractmethod
end

"""
    get_mass(op::TransientFEOperator{MassLinearODE}) -> Function

Return the mass bilinear form of the `TransientFEOperator`.
"""
function get_mass(op::TransientFEOperator{MassLinearODE})
  @abstractmethod
end

"""
    get_forms(op::TransientFEOperator{LinearODE}) -> Function

Return the bilinear forms of the `TransientFEOperator`.
"""
function get_forms(op::TransientFEOperator{LinearODE})
  @abstractmethod
end

# Broken FESpaces interface
const res_jac_on_alg_op_msg = """
For now, the residual and jacobians cannot be directly computed on a
`TransientFEOperator`. They have to be evaluated on the corresponding
algebraic operator, which is an `ODEOperator`.

This is because the `ODEOperator` works with vectors and it is optimised to
take advantage of constant jacobians.
"""

function Algebra.allocate_residual(
  op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  cache
)
  error(res_jac_on_alg_op_msg)
end

function Algebra.residual(
  op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  cache
)
  error(res_jac_on_alg_op_msg)
end

function Algebra.residual!(
  r::AbstractVector, op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  cache
)
  error(res_jac_on_alg_op_msg)
end

function Algebra.allocate_jacobian(
  op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  cache
)
  error(res_jac_on_alg_op_msg)
end

function Algebra.jacobian(
  op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  k::Integer, γ::Real,
  cache
)
  error(res_jac_on_alg_op_msg)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  k::Integer, γ::Real,
  cache
)
  error(res_jac_on_alg_op_msg)
end

function jacobians!(
  J::AbstractMatrix, op::TransientFEOperator,
  t::Real, uh::TransientCellField,
  γs::Tuple{Vararg{Real}},
  cache
)
  error(res_jac_on_alg_op_msg)
end

#############################
# TransientFEOpFromWeakForm #
#############################
"""
    struct TransientFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form.
"""
struct TransientFEOpFromWeakForm <: TransientFEOperator{NonlinearODE}
  res::Function
  jacs::Tuple{Vararg{Function}}
  jacs_constant::Tuple{Vararg{Bool}}
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

get_res(op::TransientFEOpFromWeakForm) = op.res

get_jacs(op::TransientFEOpFromWeakForm) = op.jacs

# Default constructors
function TransientFEOperator(
  res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1)
)
  if trial isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end

  function jac_0(t, u, du, v)
    function res_0(y)
      u0 = TransientCellFieldType(y, u.derivatives)
      res(t, u0, v)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for k in 1:order
    function jac_k(t, u, duk, v)
      function res_k(y)
        derivatives = (u.derivatives[1:k-1]..., y, u.derivatives[k+1:end]...)
        uk = TransientCellFieldType(u.cellfield, derivatives)
        res(t, uk, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  TransientFEOperator(res, jacs..., trial, test; jacs_constant)
end

function TransientFEOperator(
  res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  jacs_constant::NTuple{2,Bool}=ntuple(_ -> false, 2)
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientFEOpFromWeakForm(
    res, (jac, jac_t), jacs_constant,
    assembler, trial, test, 1
  )
end

function TransientFEOperator(
  res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  jacs_constant::NTuple{3,Bool}=ntuple(_ -> false, 3)
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientFEOpFromWeakForm(
    res, (jac, jac_t, jac_tt), jacs_constant,
    assembler, trial, test, 2
  )
end

#######################################
# TransientMassLinearFEOpFromWeakForm #
#######################################
"""
    struct TransientMassLinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```math
mass(t, u, v) + res(t, u, v) = 0
```, where `mass` is linear in `∂t^N(u)` (`N` is the order of the operator) and
`res` does not depend on `∂t^N(u)` (`mass` and `res` can depend on lower-order
time-derivatives of `u` by applying the `∂t` operator).
"""
struct TransientMassLinearFEOpFromWeakForm <: TransientFEOperator{MassLinearODE}
  mass::Function
  res::Function
  jacs::Tuple{Vararg{Function}}
  jacs_constant::Tuple{Vararg{Bool}}
  forcing_constant::Bool
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

get_res(op::TransientMassLinearFEOpFromWeakForm) = op.res

get_jacs(op::TransientMassLinearFEOpFromWeakForm) = op.jacs

get_mass(op::TransientMassLinearFEOpFromWeakForm) = op.mass

# Default constructors
function TransientMassLinearFEOperator(
  mass::Function, res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1),
  forcing_constant::Bool=false
)
  if trial isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end

  function jac_0(t, u, du, v)
    function res_0(y)
      u0 = TransientCellFieldType(y, u.derivatives)
      mass(t, u0, v) + res(t, u0, v)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for k in 1:order-1
    function jac_k(t, u, duk, v)
      function res_k(y)
        derivatives = (u.derivatives[1:k-1]..., y, u.derivatives[k+1:end]...)
        uk = TransientCellFieldType(u.cellfield, derivatives)
        mass(t, uk, v) + res(t, uk, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  function jac_N(t, u, duN, v)
    function res_N(y)
      derivatives = (u.derivatives[1:order-1]..., y)
      uN = TransientCellFieldType(u.cellfield, derivatives)
      mass(t, uN, v)
    end
    jacobian(res_N, u.derivatives[order])
  end
  jacs = (jacs..., jac_N)

  TransientMassLinearFEOperator(
    mass, res, jacs..., trial, test;
    jacs_constant, forcing_constant
  )
end

function TransientMassLinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  jacs_constant::NTuple{2,Bool}=ntuple(_ -> false, 2),
  forcing_constant::Bool=false
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientMassLinearFEOpFromWeakForm(
    mass, res, (jac, jac_t), jacs_constant, forcing_constant,
    assembler, trial, test, 1
  )
end

function TransientMassLinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  jacs_constant::NTuple{3,Bool}=ntuple(_ -> false, 3),
  forcing_constant::Bool=false
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientMassLinearFEOpFromWeakForm(
    mass, res, (jac, jac_t, jac_tt), jacs_constant, forcing_constant,
    assembler, trial, test, 2
  )
end

###################################
# TransientLinearFEOpFromWeakForm #
###################################
# TODO Fix notations
"""
    struct TransientLinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```math
∑_{0 ≤ k ≤ N} form_k(t, ∂t^k(u), v) + res(t, v) = 0
```, where `N` is the order of the operator, `form_k` is linear in `∂t^k(u)`
and does not depend on the other time derivatives of `u`.
"""
struct TransientLinearFEOpFromWeakForm <: TransientFEOperator{LinearODE}
  forms::Tuple{Vararg{Function}}
  res::Function
  jacs::Tuple{Vararg{Function}}
  jacs_constant::Tuple{Vararg{Bool}}
  forcing_constant::Bool
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

get_res(op::TransientLinearFEOpFromWeakForm) = (t, u, v) -> op.res(t, v)

get_jacs(op::TransientLinearFEOpFromWeakForm) = op.jacs

get_forms(op::TransientLinearFEOpFromWeakForm) = op.forms

# Default constructors
function TransientLinearFEOperator(
  forms::Tuple{Vararg{Function}}, res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1),
  forcing_constant::Bool=false
)
  if trial isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end

  function jac_0(t, u, du, v)
    function res_0(y)
      u0 = TransientCellFieldType(y, u.derivatives)
      forms[order+1](t, u0, v)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for k in 1:order
    function jac_k(t, u, duk, v)
      function res_k(y)
        derivatives = (u.derivatives[1:k-1]..., y, u.derivatives[k+1:end]...)
        uk = TransientCellFieldType(u.cellfield, derivatives)
        forms[order+1-k](t, uk, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  TransientLinearFEOperator(
    forms, res, jacs..., trial, test;
    jacs_constant, forcing_constant
  )
end

function TransientLinearFEOperator(
  mass::Function, stiffness::Function, res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  jacs_constant::NTuple{2,Bool}=ntuple(_ -> false, 2),
  forcing_constant::Bool=false
)
  TransientLinearFEOperator(
    (mass, stiffness), res,
    jac, jac_t, trial, test;
    jacs_constant, forcing_constant
  )
end

function TransientLinearFEOperator(
  forms::NTuple{2,Function}, res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  jacs_constant::NTuple{2,Bool}=ntuple(_ -> false, 2),
  forcing_constant::Bool=false
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientLinearFEOpFromWeakForm(
    forms, res, (jac, jac_t), jacs_constant, forcing_constant,
    assembler, trial, test, 1
  )
end

function TransientLinearFEOperator(
  mass::Function, stiffness::Function, damping::Function, res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  jacs_constant::NTuple{3,Bool}=ntuple(_ -> false, 3),
  forcing_constant::Bool=false
)
  TransientLinearFEOperator(
    (mass, stiffness, damping), res,
    jac, jac_t, jac_tt, trial, test;
    jacs_constant, forcing_constant
  )
end

function TransientLinearFEOperator(
  forms::NTuple{3,Function}, res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  jacs_constant::NTuple{3,Bool}=ntuple(_ -> false, 3),
  forcing_constant::Bool=false
)
  assembler = SparseMatrixAssembler(trial, test)
  TransientLinearFEOpFromWeakForm(
    forms, res, (jac, jac_t, jac_tt), jacs_constant, forcing_constant,
    assembler, trial, test, 2
  )
end

####################
# Common functions #
####################
const TransientFEOperatorTypes = Union{
  TransientFEOpFromWeakForm,
  TransientMassLinearFEOpFromWeakForm,
  TransientLinearFEOpFromWeakForm
}
const TransientAbstractMassLinearFEOperatorTypes = Union{
  TransientMassLinearFEOpFromWeakForm,
  TransientLinearFEOpFromWeakForm
}

# FEOperator interface
FESpaces.get_test(op::TransientFEOperatorTypes) = op.test

FESpaces.get_trial(op::TransientFEOperatorTypes) = op.trial

# ODEOperator interface
Polynomials.get_order(op::TransientFEOperatorTypes) = op.order

# TransientFEOperator interface
get_assembler(op::TransientFEOperatorTypes) = op.assembler
function is_jacobian_constant(op::TransientFEOperatorTypes, k::Integer)
  op.jacs_constant[k+1]
end

function is_forcing_constant(op::TransientAbstractMassLinearFEOperatorTypes)
  op.forcing_constant
end

########
# Test #
########
"""
    test_transient_fe_operator(
      op::TransientFEOperator,
      t::Real, uh::TransientCellField
    ) -> Bool

Test the interface of `TransientFEOperator` specializations.
"""
function test_transient_fe_operator(
  op::TransientFEOperator,
  t::Real, uh::TransientCellField
)
  U = get_trial(op)
  Ut = U(t)
  @test Ut isa FESpace

  V = get_test(op)
  @test V isa FESpace

  us = (get_free_dof_values(uh.cellfield),)
  for derivative in uh.derivatives
    us = (us..., get_free_dof_values(derivative))
  end

  ode_op = get_algebraic_operator(op)
  @test ode_op isa ODEOperator

  test_ode_operator(ode_op, t, us, t, us)

  true
end
