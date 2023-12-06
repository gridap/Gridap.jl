"""
    abstract type TransientFEOperator <: GridapType end

Transient version of `FEOperator` corresponding to a residual of the form
```math
residual(t, u, v) = 0,
```
where `residual` is linear in `v`. Time derivatives of `u` can be included by
using the `∂t` operator.

# Important
For now, the residual and jacobians cannot be directly computed on a
`TransientFEOperator`. They have to be evaluated on the corresponding algebraic
operator, which is an `ODEOperator`. As such, `TransientFEOperator` is not
exactly a subtype of `FEOperator`, but rather at the intersection of
`FEOperator` and `ODEOperator`.

# Mandatory
- [`get_test(feop)`](@ref)
- [`get_trial(feop)`](@ref)
- [`get_order(feop)`](@ref)
- [`get_assembler(feop)`](@ref)
- [`get_res(feop::TransientFEOperator)`](@ref)
- [`get_jacs(feop::TransientFEOperator)`](@ref)
- [`get_mass(feop::TransientFEOperator{<:AbstractQuasilinearODE})`](@ref)
- [`get_forms(feop::TransientFEOperator{<:AbstractLinearODE})`](@ref)

# Optional
- [`allocate_feopcache(feop)`](@ref)
- [`update_feopcache!(feopcache, feop, t)`](@ref)
- [`get_algebraic_operator(feop)`](@ref)
- [`is_jacobian_constant(feop, k)`](@ref)
- [`is_residual_constant(feop)`](@ref)
"""
abstract type TransientFEOperator{C<:ODEOperatorType} <: GridapType end

"""
    ODEOperatorType(::Type{<:TransientFEOperator}) -> ODEOperatorType

Return the `ODEOperatorType` of the `TransientFEOperator`.
"""
ODEOperatorType(feop::TransientFEOperator) = ODEOperatorType(typeof(feop))
ODEOperatorType(::Type{<:TransientFEOperator{C}}) where {C} = C

# FEOperator interface
function FESpaces.get_test(feop::TransientFEOperator)
  @abstractmethod
end

function FESpaces.get_trial(feop::TransientFEOperator)
  @abstractmethod
end

function FESpaces.get_algebraic_operator(feop::TransientFEOperator)
  ODEOpFromFEOp(feop)
end

# ODEOperator interface
function Polynomials.get_order(feop::TransientFEOperator)
  @abstractmethod
end

function is_jacobian_constant(feop::TransientFEOperator, k::Integer)
  false
end

function is_residual_constant(feop::TransientFEOperator)
  false
end

# TransientFEOperator interface
"""
    allocate_feopcache(feop::TransientFEOperator)

Allocate the cache of the `TransientFEOperator`.
"""
function allocate_feopcache(feop::TransientFEOperator)
  nothing
end

"""
    update_feopcache!(feopcache, feop::TransientFEOperator, t::Real)

Update the cache of the `TransientFEOperator` at time `t`.
"""
function update_feopcache!(feopcache, feop::TransientFEOperator, t::Real)
  feopcache
end

"""
    get_assembler(feop::TransientFEOperator) -> Assembler

Return the assembler of the `TransientFEOperator`.
"""
function get_assembler(feop::TransientFEOperator)
  @abstractmethod
end

"""
    get_res(feop::TransientFEOperator) -> Function

Return the lowest-order element in the decomposition of the residual of the
`ODEOperator`:
* In the general case, return the whole residual,
* For an `AbstractQuasilinearODE`, return the residual excluding the mass term,
* For an `AbstractLinearODE`, return the forcing term.
"""
function get_res(feop::TransientFEOperator)
  @abstractmethod
end

"""
    get_jacs(feop::TransientFEOperator) -> Tuple{Vararg{Function}}

Return the jacobians of the `TransientFEOperator`.
"""
function get_jacs(feop::TransientFEOperator)
  @abstractmethod
end

"""
    get_mass(feop::TransientFEOperator{<:AbstractQuasilinearODE}) -> Function

Return the mass bilinear form of the `TransientFEOperator`.
"""
function get_mass(feop::TransientFEOperator{<:AbstractQuasilinearODE})
  @abstractmethod
end

function get_mass(feop::TransientFEOperator{<:AbstractLinearODE})
  get_forms(feop)[end]
end

"""
    get_forms(feop::TransientFEOperator{<:AbstractLinearODE}) -> Function

Return the bilinear forms of the `TransientFEOperator`.
"""
function get_forms(feop::TransientFEOperator{<:AbstractLinearODE})
  @abstractmethod
end

# Broken FESpaces interface
const res_jac_on_transient_feop_msg = """
For now, the residual and jacobians cannot be directly computed on a
`TransientFEOperator`. They have to be evaluated on the corresponding
algebraic operator, which is an `ODEOperator`.

This is because the `ODEOperator` works with vectors and it is optimised to
take advantage of constant jacobians.
"""

function Algebra.allocate_residual(
  feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

function Algebra.residual(
  feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

function Algebra.residual!(
  r::AbstractVector, feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

function Algebra.allocate_jacobian(
  feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

function Algebra.jacobian(
  feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  k::Integer, γ::Real,
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

function Algebra.jacobian!(
  J::AbstractMatrix, feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  k::Integer, γ::Real,
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

function jacobians!(
  J::AbstractMatrix, feop::TransientFEOperator,
  t::Real, uh::TransientCellField,
  γs::Tuple{Vararg{Real}},
  feopcache
)
  error(res_jac_on_transient_feop_msg)
end

###########################
# TransientIMEXFEOperator #
###########################
"""
    struct TransientIMEXFEOperator <: TransientFEOperator end

Implicit-Explicit version of `TransientFEOperator` whose residual can be
decomposed into
```math
residual(t, u, v) = implicit_residual(t, u, v)
                  + explicit_residual(t, u, v),
```
where `implicit_residual` and `explicit_residual` are linear in `v` and the
explicit residual has one order less than the implicit residual.
"""
struct TransientIMEXFEOperator{Cim,Cex} <: TransientFEOperator{Cim}
  im_feop::TransientFEOperator{Cim}
  ex_feop::TransientFEOperator{Cex}

  function TransientIMEXFEOperator(
    im_feop::TransientFEOperator,
    ex_feop::TransientFEOperator
  )
    msg_spaces = """
    The implicit and explicit `TransientFEOperator` of a
    `TransientIMEXFEOperator` must be defined on the same trial and test spaces.
    """
    @assert get_trial(im_feop) == get_trial(ex_feop) msg_spaces
    @assert get_test(im_feop) == get_test(ex_feop) msg_spaces

    msg_order = """
    The explicit `TransientFEOperator` of a `TransientIMEXFEOperator` must have
    one order less than the implicit `TransientFEOperator`.
    """
    @assert get_order(im_feop) == get_order(ex_feop) + 1 msg_order
    Cim = ODEOperatorType(im_feop)
    Cex = ODEOperatorType(ex_feop)
    new{Cim,Cex}(im_feop, ex_feop)
  end
end

# TransientFEOperator interface
# Only these function need to be implemented because all other functions of the
# interface are going to be called on the implicit and explicit
# `ODEOpFromFEOp`s within the `IMEXODEOperator` interface, and in turn called
# on the implicit and explicit `TransientFEOperator`s separately
function FESpaces.get_trial(feop::TransientIMEXFEOperator)
  feop.im_feop.trial
end

function FESpaces.get_test(feop::TransientIMEXFEOperator)
  feop.im_feop.tes
end

function FESpaces.get_algebraic_operator(feop::TransientIMEXFEOperator)
  im_odeop = ODEOpFromFEOp(feop.im_feop)
  ex_odeop = ODEOpFromFEOp(feop.ex_feop)
  GenericIMEXODEOperator(im_odeop, ex_odeop)
end

#############################
# TransientFEOpFromWeakForm #
#############################
"""
    struct TransientFEOpFromWeakForm <: TransientFEOperator end

Generic `TransientFEOperator` constructed from the weak formulation of a
partial differential equation.
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

# Constructors
function TransientFEOperator(
  res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1)
)
  function jac_0(t, u, du, v)
    function res_0(y)
      u0 = TransientCellField(y, u.derivatives)
      res(t, u0, v)
    end
    jacobian(res_0, u.cellfield)
  end
  jacs = (jac_0,)

  for k in 1:order
    function jac_k(t, u, duk, v)
      function res_k(y)
        derivatives = (u.derivatives[1:k-1]..., y, u.derivatives[k+1:end]...)
        uk = TransientCellField(u.cellfield, derivatives)
        res(t, uk, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  TransientFEOperator(res, jacs, trial, test; jacs_constant)
end

function TransientFEOperator(
  res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, length(jacs))
)
  order = length(jacs) - 1
  assembler = SparseMatrixAssembler(trial, test)
  TransientFEOpFromWeakForm(
    res, jacs, jacs_constant,
    assembler, trial, test, order
  )
end

# TransientFEOperator interface (rest in # Common functions)
get_res(feop::TransientFEOpFromWeakForm) = feop.res

get_jacs(feop::TransientFEOpFromWeakForm) = feop.jacs

########################################
# TransientQuasilinearFEOpFromWeakForm #
########################################
"""
    struct TransientQuasilinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```math
residual(t, u, v) = mass(t, u, ∂t^N[u] v) + res(t, u, v) = 0.
```
Let `N` be the order of the operator. We impose the following conditions:
* `mass` is linear in the `N`-th-order time derivative of `u`,
* `res` has order `N-1`,
* both `mass` and `res` are linear in `v`.

For convenience, the mass matrix has to be specified as a function of `u` for
the nonlinear part, and `∂t^N[u]`.
"""
struct TransientQuasilinearFEOpFromWeakForm <: TransientFEOperator{QuasilinearODE}
  mass::Function
  res::Function
  jacs::Tuple{Vararg{Function}}
  jacs_constant::Tuple{Vararg{Bool}}
  residual_constant::Bool
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructors
function TransientQuasilinearFEOperator(
  mass::Function, res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1),
  residual_constant::Bool=false
)
  jacs = ()
  if order > 0
    function jac_0(t, u, du, v)
      function res_0(y)
        u0 = TransientCellField(y, u.derivatives)
        ∂tNu0 = u0
        for _ in 1:order
          ∂tNu0 = ∂t(∂tNu0)
        end
        mass(t, u0, ∂tNu0, v) + res(t, u0, v)
      end
      jacobian(res_0, u.cellfield)
    end
    jacs = (jacs..., jac_0)
  end

  for k in 1:order-1
    function jac_k(t, u, duk, v)
      function res_k(y)
        derivatives = (u.derivatives[1:k-1]..., y, u.derivatives[k+1:end]...)
        uk = TransientCellField(u.cellfield, derivatives)
        ∂tNuk = uk
        for _ in 1:order
          ∂tNuk = ∂t(∂tNuk)
        end
        mass(t, uk, ∂tNuk, v) + res(t, uk, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  function jac_N(t, u, duN, v)
    function res_N(y)
      derivatives = (u.derivatives[1:end-1]..., y)
      uN = TransientCellField(u.cellfield, derivatives)
      mass(t, uN, y, v)
    end
    jacobian(res_N, u.derivatives[end])
  end
  jacs = (jacs..., jac_N)

  TransientQuasilinearFEOperator(
    mass, res, jacs, trial, test;
    jacs_constant, residual_constant
  )
end

function TransientQuasilinearFEOperator(
  mass::Function, res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, length(jacs)),
  residual_constant::Bool=false
)
  order = length(jacs) - 1
  assembler = SparseMatrixAssembler(trial, test)
  TransientQuasilinearFEOpFromWeakForm(
    mass, res, jacs, jacs_constant, residual_constant,
    assembler, trial, test, order
  )
end

# TransientFEOperator interface (rest in # Common functions)
get_res(feop::TransientQuasilinearFEOpFromWeakForm) = feop.res

get_jacs(feop::TransientQuasilinearFEOpFromWeakForm) = feop.jacs

get_mass(feop::TransientQuasilinearFEOpFromWeakForm) = feop.mass

########################################
# TransientQuasilinearFEOpFromWeakForm #
########################################
"""
    struct TransientSemilinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```math
residual(t, u, v) = mass(t, ∂t^N[u], v) + res(t, u, v) = 0.
```
Let `N` be the order of the operator. We impose the following conditions:
* `mass` is linear in the `N`-th-order time derivative of `u`,
* `res` has order `N-1`,
* both `mass` and `res` are linear in `v`.

For convenience, the mass matrix has to be specified as a function of
`∂t^N[u]`, i.e. as a linear form.
"""
struct TransientSemilinearFEOpFromWeakForm <: TransientFEOperator{SemilinearODE}
  mass::Function
  res::Function
  jacs::Tuple{Vararg{Function}}
  jacs_constant::Tuple{Vararg{Bool}}
  residual_constant::Bool
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructors
function TransientSemilinearFEOperator(
  mass::Function, res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1),
  residual_constant::Bool=false
)
  # When the operator is semilinear, the mass term can be omitted in the
  # computation of the other jacobians.
  jacs = ()
  if order > 0
    function jac_0(t, u, du, v)
      function res_0(y)
        u0 = TransientCellField(y, u.derivatives)
        res(t, u0, v)
      end
      jacobian(res_0, u.cellfield)
    end
    jacs = (jacs..., jac_0)
  end

  for k in 1:order-1
    function jac_k(t, u, duk, v)
      function res_k(y)
        derivatives = (u.derivatives[1:k-1]..., y, u.derivatives[k+1:end]...)
        uk = TransientCellField(u.cellfield, derivatives)
        res(t, uk, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  # When the operator is semilinear, the jacobian of the residual w.r.t. the
  # highest-order term is simply the mass term.
  jac_N(t, u, duN, v) = mass(t, duN, v)
  jacs = (jacs..., jac_N)

  TransientSemilinearFEOperator(
    mass, res, jacs, trial, test;
    jacs_constant, residual_constant
  )
end

function TransientSemilinearFEOperator(
  mass::Function, res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, length(jacs)),
  residual_constant::Bool=false,
)
  order = length(jacs) - 1
  assembler = SparseMatrixAssembler(trial, test)
  TransientSemilinearFEOpFromWeakForm(
    mass, res, jacs, jacs_constant, residual_constant,
    assembler, trial, test, order
  )
end

# TransientFEOperator interface (rest in # Common functions)
get_res(feop::TransientSemilinearFEOpFromWeakForm) = feop.res

get_jacs(feop::TransientSemilinearFEOpFromWeakForm) = feop.jacs

get_mass(feop::TransientSemilinearFEOpFromWeakForm) = feop.mass

###################################
# TransientLinearFEOpFromWeakForm #
###################################
"""
    struct TransientLinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```math
residual(t, u, v) = ∑_{0 ≤ k ≤ N} form_k(t, ∂t^k[u], v) + res(t, v) = 0,
```
where `N` is the order of the operator, `form_k` is linear in `∂t^k[u]` and
does not depend on the other time derivatives of `u`, and the `form_k` and
`res` are linear in `v`.

For convenience, the form corresponding to order `k` has to be written as a
function of `∂t^k[u]`, i.e. as a linear form, and the residual as a function
of `t` and `v` only.
"""
struct TransientLinearFEOpFromWeakForm <: TransientFEOperator{LinearODE}
  forms::Tuple{Vararg{Function}}
  res::Function
  jacs::Tuple{Vararg{Function}}
  jacs_constant::Tuple{Vararg{Bool}}
  residual_constant::Bool
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructors
function TransientLinearFEOperator(
  forms::Tuple{Vararg{Function}}, res::Function,
  trial, test;
  order::Integer=1,
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, order + 1),
  residual_constant::Bool=false,
)
  # When the operator is linear, the jacobians are the forms themselves
  jacs = ntuple(k -> ((t, u, duk, v) -> forms[k](t, duk, v)), order + 1)

  TransientLinearFEOperator(
    forms, res, jacs, trial, test;
    jacs_constant, residual_constant
  )
end

function TransientLinearFEOperator(
  forms::Tuple{Vararg{Function}}, res::Function,
  jacs::Tuple{Vararg{Function}},
  trial, test;
  jacs_constant::Tuple{Vararg{Bool}}=ntuple(_ -> false, length(jacs)),
  residual_constant::Bool=false
)
  order = length(jacs) - 1
  assembler = SparseMatrixAssembler(trial, test)
  TransientLinearFEOpFromWeakForm(
    forms, res, jacs, jacs_constant, residual_constant,
    assembler, trial, test, order
  )
end

# TransientFEOperator interface (rest in # Common functions)
get_res(feop::TransientLinearFEOpFromWeakForm) = (t, u, v) -> feop.res(t, v)

get_jacs(feop::TransientLinearFEOpFromWeakForm) = feop.jacs

get_forms(feop::TransientLinearFEOpFromWeakForm) = feop.forms

####################
# Common functions #
####################
const TransientFEOperatorTypes = Union{
  TransientFEOpFromWeakForm,
  TransientQuasilinearFEOpFromWeakForm,
  TransientSemilinearFEOpFromWeakForm,
  TransientLinearFEOpFromWeakForm
}
const TransientAbstractQuasilinearFEOperatorTypes = Union{
  TransientQuasilinearFEOpFromWeakForm,
  TransientSemilinearFEOpFromWeakForm,
  TransientLinearFEOpFromWeakForm
}

# FEOperator interface
FESpaces.get_test(feop::TransientFEOperatorTypes) = feop.test

FESpaces.get_trial(feop::TransientFEOperatorTypes) = feop.trial

# ODEOperator interface
Polynomials.get_order(feop::TransientFEOperatorTypes) = feop.order

# TransientFEOperator interface
get_assembler(feop::TransientFEOperatorTypes) = feop.assembler

function is_jacobian_constant(feop::TransientFEOperatorTypes, k::Integer)
  feop.jacs_constant[k+1]
end

function is_residual_constant(feop::TransientAbstractQuasilinearFEOperatorTypes)
  feop.residual_constant
end

########
# Test #
########
"""
    test_transient_fe_operator(
      feop::TransientFEOperator,
      t::Real, uh::TransientCellField
    ) -> Bool

Test the interface of `TransientFEOperator` specializations.
"""
function test_transient_fe_operator(
  feop::TransientFEOperator,
  t::Real, uh::TransientCellField
)
  U = get_trial(feop)
  Ut = U(t)
  @test Ut isa FESpace

  V = get_test(feop)
  @test V isa FESpace

  us = (get_free_dof_values(uh.cellfield),)
  for derivative in uh.derivatives
    us = (us..., get_free_dof_values(derivative))
  end

  odeop = get_algebraic_operator(feop)
  @test odeop isa ODEOperator

  test_ode_operator(odeop, t, us)

  true
end
