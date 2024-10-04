"""
    abstract type TransientFEOperator <: GridapType end

Transient version of `FEOperator` corresponding to a residual of the form
```
residual(t, u, v) = 0,
```
where `residual` is linear in `v`. Time derivatives of `u` can be included by
using the `∂t` operator.

# Important
For now, the residual and jacobians cannot be directly computed on a
`TransientFEOperator`. They have to be evaluated on the corresponding
algebraic operator, which is an `ODEOperator`. As such, `TransientFEOperator`
is not exactly a subtype of `FEOperator`, but rather at the intersection of
`FEOperator` and `ODEOperator`. This is because the `ODEOperator` works with
vectors and it is optimised to take advantage of constant forms.

# Mandatory
- [`get_test(tfeop)`](@ref)
- [`get_trial(tfeop)`](@ref)
- [`get_order(tfeop)`](@ref)
- [`get_res(tfeop::TransientFEOperator)`](@ref)
- [`get_jacs(tfeop::TransientFEOperator)`](@ref)
- [`get_forms(tfeop::TransientFEOperator)`](@ref)
- [`get_assembler(tfeop)`](@ref)

# Optional
- [`get_algebraic_operator(tfeop)`](@ref)
- [`get_num_forms(tfeop::TransientFEOperator)`](@ref)
- [`is_form_constant(tfeop, k)`](@ref)
- [`allocate_tfeopcache(tfeop)`](@ref)
- [`update_tfeopcache!(tfeopcache, tfeop, t)`](@ref)
"""
abstract type TransientFEOperator{T<:ODEOperatorType} <: GridapType end

"""
    ODEOperatorType(::Type{<:TransientFEOperator}) -> ODEOperatorType

Return the `ODEOperatorType` of the `TransientFEOperator`.
"""
ODEOperatorType(::TransientFEOperator{T}) where {T} = T
ODEOperatorType(::Type{<:TransientFEOperator{T}}) where {T} = T

# FEOperator interface
function FESpaces.get_test(tfeop::TransientFEOperator)
  @abstractmethod
end

function FESpaces.get_trial(tfeop::TransientFEOperator)
  @abstractmethod
end

function FESpaces.get_algebraic_operator(tfeop::TransientFEOperator)
  ODEOpFromTFEOp(tfeop)
end

# ODEOperator interface
function Polynomials.get_order(tfeop::TransientFEOperator)
  @abstractmethod
end

# TransientFEOperator interface
"""
    get_res(tfeop::TransientFEOperator) -> Function

Return the lowest-order element in the decomposition of the residual of the
`ODEOperator`:
* In the general case, return the whole residual,
* For an `AbstractQuasilinearODE`, return the residual excluding the mass term,
* For an `AbstractLinearODE`, return the forcing term.
"""
function get_res(tfeop::TransientFEOperator)
  @abstractmethod
end

"""
    get_jacs(tfeop::TransientFEOperator) -> Tuple{Vararg{Function}}

Return the jacobians of the `TransientFEOperator`.
"""
function get_jacs(tfeop::TransientFEOperator)
  @abstractmethod
end

"""
    get_num_forms(tfeop::TransientFEOperator) -> Integer

Return the number of bilinear forms of the `TransientFEOperator`. See
[`get_forms`](@ref).
"""
function get_num_forms(tfeop::TransientFEOperator)
  0
end

function get_num_forms(tfeop::TransientFEOperator{<:AbstractQuasilinearODE})
  1
end

function get_num_forms(tfeop::TransientFEOperator{<:AbstractLinearODE})
  get_order(tfeop) + 1
end

"""
    get_forms(tfeop::TransientFEOperator) -> Function

Return the bilinear forms of the `TransientFEOperator`:
* For a general transient FE operator, return nothing,
* For a quasilinear transient FE operator, return the mass matrix,
* For a linear transient FE operator, return all the linear forms.
"""
function get_forms(tfeop::TransientFEOperator)
  ()
end

function get_forms(tfeop::TransientFEOperator{<:AbstractQuasilinearODE})
  @abstractmethod
end

"""
    is_form_constant(tfeop::TransientFEOperator, k::Integer) -> Bool

Indicate whether the bilinear form of the `TransientFEOperator` corresponding
to the `k`-th-order time derivative of `u` is constant with respect to `t`.
"""
function is_form_constant(tfeop::TransientFEOperator, k::Integer)
  false
end

"""
    get_assembler(tfeop::TransientFEOperator) -> Assembler

Return the assembler of the `TransientFEOperator`.
"""
function get_assembler(tfeop::TransientFEOperator)
  @abstractmethod
end

"""
    allocate_tfeopcache(
      tfeop::TransientFEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}
    ) -> CacheType

Allocate the cache of the `TransientFEOperator`.
"""
function allocate_tfeopcache(
  tfeop::TransientFEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}
)
  nothing
end

"""
    update_tfeopcache!(tfeopcache, tfeop::TransientFEOperator, t::Real) -> CacheType

Update the cache of the `TransientFEOperator` at time `t`.
"""
function update_tfeopcache!(tfeopcache, tfeop::TransientFEOperator, t::Real)
  tfeopcache
end

# Broken FESpaces interface
const res_jac_on_transient_tfeop_msg = """
For now, the residual and jacobians cannot be directly computed on a
`TransientFEOperator`. They have to be evaluated on the corresponding
algebraic operator, which is an `ODEOperator`.

This is because the `ODEOperator` works with vectors and it is optimised to
take advantage of constant jacobians.
"""

function Algebra.allocate_residual(tfeop::TransientFEOperator, u)
  @unreachable res_jac_on_transient_tfeop_msg
end

function Algebra.residual!(r::AbstractVector, tfeop::TransientFEOperator, u)
  @unreachable res_jac_on_transient_tfeop_msg
end

function Algebra.residual(tfeop::TransientFEOperator, u)
  @unreachable res_jac_on_transient_tfeop_msg
end

function Algebra.allocate_jacobian(tfeop::TransientFEOperator, u)
  @unreachable res_jac_on_transient_tfeop_msg
end

function Algebra.jacobian!(J::AbstractMatrix, tfeop::TransientFEOperator, u)
  @unreachable res_jac_on_transient_tfeop_msg
end

function Algebra.jacobian(tfeop::TransientFEOperator, u)
  @unreachable res_jac_on_transient_tfeop_msg
end

const default_linear_msg = """
For an operator of order zero, the definitions of quasilinear, semilinear and
linear coincide. Defaulting to linear.
"""

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
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructor with manual jacobians
function TransientFEOperator(
  res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  order = length(jacs) - 1
  TransientFEOpFromWeakForm(
    res, jacs,
    assembler, trial, test, order
  )
end

# Constructors with flat arguments (orders 0, 1, 2)
function TransientFEOperator(
  res::Function,
  jac::Function,
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientFEOperator(
    res, (jac,), trial, test;
    assembler
  )
end

function TransientFEOperator(
  res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientFEOperator(
    res, (jac, jac_t), trial, test;
    assembler
  )
end

function TransientFEOperator(
  res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientFEOperator(
    res, (jac, jac_t, jac_tt), trial, test;
    assembler
  )
end

# Constructor with automatic jacobians
function TransientFEOperator(
  res::Function,
  trial, test;
  order::Integer=1,
  assembler=SparseMatrixAssembler(trial, test)
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

  TransientFEOperator(
    res, jacs, trial, test;
    assembler
  )
end

# TransientFEOperator interface
FESpaces.get_test(tfeop::TransientFEOpFromWeakForm) = tfeop.test

FESpaces.get_trial(tfeop::TransientFEOpFromWeakForm) = tfeop.trial

Polynomials.get_order(tfeop::TransientFEOpFromWeakForm) = tfeop.order

get_res(tfeop::TransientFEOpFromWeakForm) = tfeop.res

get_jacs(tfeop::TransientFEOpFromWeakForm) = tfeop.jacs

get_assembler(tfeop::TransientFEOpFromWeakForm) = tfeop.assembler

########################################
# TransientQuasilinearFEOpFromWeakForm #
########################################
"""
    struct TransientQuasilinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```
residual(t, u, v) = mass(t, u, ∂t^N[u], v) + res(t, u, v) = 0.
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
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructor with manual jacobians
function TransientQuasilinearFEOperator(
  mass::Function, res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  order = length(jacs) - 1
  if order == 0
    @warn default_linear_msg
    return TransientLinearFEOperator(
      (mass,), res, jacs, trial, test;
      assembler
    )
  end

  TransientQuasilinearFEOpFromWeakForm(
    mass, res, jacs,
    assembler, trial, test, order
  )
end

# Constructor with flat arguments (orders 0, 1, 2)
function TransientQuasilinearFEOperator(
  mass::Function, res::Function,
  jac::Function,
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  @warn default_linear_msg
  TransientLinearFEOperator(
    mass, res, jac, trial, test;
    assembler
  )
end

function TransientQuasilinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientQuasilinearFEOperator(
    mass, res, (jac, jac_t), trial, test;
    assembler
  )
end

function TransientQuasilinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientQuasilinearFEOperator(
    mass, res, (jac, jac_t, jac_tt), trial, test;
    assembler
  )
end

# Constructor with automatic jacobians
function TransientQuasilinearFEOperator(
  mass::Function, res::Function,
  trial, test;
  order::Integer=1,
  assembler=SparseMatrixAssembler(trial, test)
)
  if order == 0
    @warn default_linear_msg
    return TransientLinearFEOperator(
      mass, res, trial, test;
      assembler
    )
  end

  jacs = ()
  if order > 0
    function jac_0(t, u, du, v)
      function res_0(y)
        u0 = TransientCellField(y, u.derivatives)
        ∂tNu0 = ∂t(u0, Val(order))
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
        u0 = TransientCellField(u.cellfield, derivatives)
        ∂tNu0 = ∂t(u0, Val(order))
        mass(t, u0, ∂tNu0, v) + res(t, u0, v)
      end
      jacobian(res_k, u.derivatives[k])
    end
    jacs = (jacs..., jac_k)
  end

  # When the operator is quasilinear, the jacobian of the residual w.r.t. the
  # highest-order term is simply the mass term.
  jac_N(t, u, duN, v) = mass(t, u, duN, v)
  jacs = (jacs..., jac_N)

  TransientQuasilinearFEOperator(
    mass, res, jacs, trial, test;
    assembler
  )
end

# TransientFEOperator interface
FESpaces.get_test(tfeop::TransientQuasilinearFEOpFromWeakForm) = tfeop.test

FESpaces.get_trial(tfeop::TransientQuasilinearFEOpFromWeakForm) = tfeop.trial

Polynomials.get_order(tfeop::TransientQuasilinearFEOpFromWeakForm) = tfeop.order

get_res(tfeop::TransientQuasilinearFEOpFromWeakForm) = tfeop.res

get_jacs(tfeop::TransientQuasilinearFEOpFromWeakForm) = tfeop.jacs

get_forms(tfeop::TransientQuasilinearFEOpFromWeakForm) = (tfeop.mass,)

get_assembler(tfeop::TransientQuasilinearFEOpFromWeakForm) = tfeop.assembler

########################################
# TransientQuasilinearFEOpFromWeakForm #
########################################
"""
    struct TransientSemilinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```
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
  constant_mass::Bool
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructor with manual jacobians
function TransientSemilinearFEOperator(
  mass::Function, res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  constant_mass::Bool=false,
  assembler=SparseMatrixAssembler(trial, test)
)
  order = length(jacs) - 1
  if order == 0
    @warn default_linear_msg
    constant_forms = (constant_mass,)
    return TransientLinearFEOperator(
      (mass,), res, jacs, trial, test;
      constant_forms, assembler
    )
  end

  TransientSemilinearFEOpFromWeakForm(
    mass, res, jacs, constant_mass,
    assembler, trial, test, order
  )
end

# Constructor with flat arguments (orders 0, 1, 2)
function TransientSemilinearFEOperator(
  mass::Function, res::Function,
  jac::Function,
  trial, test;
  constant_mass::Bool=false,
  assembler=SparseMatrixAssembler(trial, test)
)
  @warn default_linear_msg
  constant_forms = (constant_mass,)
  TransientLinearFEOperator(
    mass, res, jac, trial, test;
    constant_forms, assembler
  )
end

function TransientSemilinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function,
  trial, test;
  constant_mass::Bool=false,
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientSemilinearFEOperator(
    mass, res, (jac, jac_t), trial, test;
    constant_mass, assembler
  )
end

function TransientSemilinearFEOperator(
  mass::Function, res::Function,
  jac::Function, jac_t::Function, jac_tt::Function,
  trial, test;
  constant_mass::Bool=false,
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientSemilinearFEOperator(
    mass, res, (jac, jac_t, jac_tt), trial, test;
    constant_mass, assembler
  )
end

# Constructor with automatic jacobians
function TransientSemilinearFEOperator(
  mass::Function, res::Function,
  trial, test;
  order::Integer=1,
  constant_mass::Bool=false,
  assembler=SparseMatrixAssembler(trial, test)
)
  if order == 0
    @warn default_linear_msg
    constant_forms = (constant_mass,)
    return TransientLinearFEOperator(
      mass, res, trial, test;
      constant_forms, assembler
    )
  end

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
    constant_mass, assembler
  )
end

# TransientFEOperator interface
FESpaces.get_test(tfeop::TransientSemilinearFEOpFromWeakForm) = tfeop.test

FESpaces.get_trial(tfeop::TransientSemilinearFEOpFromWeakForm) = tfeop.trial

Polynomials.get_order(tfeop::TransientSemilinearFEOpFromWeakForm) = tfeop.order

get_res(tfeop::TransientSemilinearFEOpFromWeakForm) = tfeop.res

get_jacs(tfeop::TransientSemilinearFEOpFromWeakForm) = tfeop.jacs

get_forms(tfeop::TransientSemilinearFEOpFromWeakForm) = (tfeop.mass,)

function is_form_constant(tfeop::TransientSemilinearFEOpFromWeakForm, k::Integer)
  (k == get_order(tfeop)) && tfeop.constant_mass
end

get_assembler(tfeop::TransientSemilinearFEOpFromWeakForm) = tfeop.assembler

###################################
# TransientLinearFEOpFromWeakForm #
###################################
"""
    struct TransientLinearFEOpFromWeakForm <: TransientFEOperator end

Transient `FEOperator` defined by a transient weak form
```
residual(t, u, v) = ∑_{0 ≤ k ≤ N} form_k(t, ∂t^k[u], v) - res(t, v) = 0,
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
  constant_forms::Tuple{Vararg{Bool}}
  assembler::Assembler
  trial::FESpace
  test::FESpace
  order::Integer
end

# Constructor with manual jacobians
function TransientLinearFEOperator(
  forms::Tuple{Vararg{Function}}, res::Function, jacs::Tuple{Vararg{Function}},
  trial, test;
  constant_forms::Tuple{Vararg{Bool}}=ntuple(_ -> false, length(forms)),
  assembler=SparseMatrixAssembler(trial, test)
)
  order = length(jacs) - 1
  TransientLinearFEOpFromWeakForm(
    forms, res, jacs, constant_forms,
    assembler, trial, test, order
  )
end

# No constructor with flat arguments: would clash with the constructors
# below with flat forms and automatic jacobians, which are more useful

# Constructor with automatic jacobians
function TransientLinearFEOperator(
  forms::Tuple{Vararg{Function}}, res::Function,
  trial, test;
  constant_forms::Tuple{Vararg{Bool}}=ntuple(_ -> false, length(forms)),
  assembler=SparseMatrixAssembler(trial, test)
)
  # When the operator is linear, the jacobians are the forms themselves
  order = length(forms) - 1
  jacs = ntuple(k -> ((t, u, duk, v) -> forms[k](t, duk, v)), order + 1)

  TransientLinearFEOperator(
    forms, res, jacs, trial, test;
    constant_forms, assembler
  )
end

# Constructor with flat forms and automatic jacobians (orders 0, 1, 2)
function TransientLinearFEOperator(
  mass::Function, res::Function,
  trial, test;
  constant_forms::NTuple{1,Bool}=(false,),
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientLinearFEOperator(
    (mass,), res, trial, test;
    constant_forms, assembler
  )
end

function TransientLinearFEOperator(
  stiffness::Function, mass::Function, res::Function,
  trial, test;
  constant_forms::NTuple{2,Bool}=(false, false),
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientLinearFEOperator(
    (stiffness, mass), res, trial, test;
    constant_forms, assembler
  )
end

function TransientLinearFEOperator(
  stiffness::Function, damping::Function, mass::Function, res::Function,
  trial, test;
  constant_forms::NTuple{3,Bool}=(false, false, false),
  assembler=SparseMatrixAssembler(trial, test)
)
  TransientLinearFEOpFromWeakForm(
    (stiffness, damping, mass), res, trial, test;
    constant_forms, assembler
  )
end

# TransientFEOperator interface
FESpaces.get_test(tfeop::TransientLinearFEOpFromWeakForm) = tfeop.test

FESpaces.get_trial(tfeop::TransientLinearFEOpFromWeakForm) = tfeop.trial

Polynomials.get_order(tfeop::TransientLinearFEOpFromWeakForm) = tfeop.order

get_res(tfeop::TransientLinearFEOpFromWeakForm) = (t, u, v) -> tfeop.res(t, v)

get_jacs(tfeop::TransientLinearFEOpFromWeakForm) = tfeop.jacs

get_forms(tfeop::TransientLinearFEOpFromWeakForm) = tfeop.forms

function is_form_constant(tfeop::TransientLinearFEOpFromWeakForm, k::Integer)
  tfeop.constant_forms[k+1]
end

get_assembler(tfeop::TransientLinearFEOpFromWeakForm) = tfeop.assembler

###########################
# TransientIMEXFEOperator #
###########################
"""
    abstract type TransientIMEXFEOperator <: TransientFEOperator end

Implicit-Explicit decomposition of a residual defining a `TransientFEOperator`:
```
residual(t, u, v) = implicit_residual(t, u, v)
                  + explicit_residual(t, u, v),
```
where
* The implicit operator defined by the implicit residual is considered stiff and is meant to be solved implicitly,
* The explicit operator defined by the explicit residual is considered non-stiff and is meant to be solved explicitly.
* Both the implicit and explicit residuals are linear in `v`.

# Important
The explicit operator must have one order less than the implicit operator, so
that the mass term of the global operator is fully contained in the implicit
operator.

# Mandatory
- [`get_imex_operators(tfeop)`](@ref)

# Optional
- [`get_test(tfeop)`](@ref)
- [`get_trial(tfeop)`](@ref)
- [`get_algebraic_operator(tfeop)`](@ref)
"""
abstract type TransientIMEXFEOperator{T<:ODEOperatorType} <: TransientFEOperator{T} end

"""
    get_imex_operators(tfeop::TransientIMEXFEOperator) -> (TransientFEOperator, TransientFEOperator)

Return the implicit and explicit parts of the `TransientIMEXFEOperator`.
"""
function get_imex_operators(tfeop::TransientIMEXFEOperator)
  @abstractmethod
end

# TransientFEOperator interface
# Only these function need to be implemented because all other functions of the
# interface are going to be called on the implicit and explicit
# `ODEOpFromFEOp`s within the `IMEXODEOperator` interface, and in turn called
# on the implicit and explicit `TransientFEOperator`s separately
function FESpaces.get_test(tfeop::TransientIMEXFEOperator)
  im_tfeop, _ = get_imex_operators(tfeop)
  get_test(im_tfeop)
end

function FESpaces.get_trial(tfeop::TransientIMEXFEOperator)
  im_tfeop, _ = get_imex_operators(tfeop)
  get_trial(im_tfeop)
end

function FESpaces.get_algebraic_operator(tfeop::TransientIMEXFEOperator)
  im_tfeop, ex_tfeop = get_imex_operators(tfeop)
  im_odeop, ex_odeop = ODEOpFromTFEOp(im_tfeop), ODEOpFromTFEOp(ex_tfeop)
  GenericIMEXODEOperator(im_odeop, ex_odeop)
end

function FESpaces.get_order(tfeop::TransientIMEXFEOperator)
  im_tfeop, _ = get_imex_operators(tfeop)
  get_order(im_tfeop)
end

# IMEX Helpers
function check_imex_compatibility(
  im_tfeop::TransientFEOperator, ex_tfeop::TransientFEOperator
)
  msg = """
  The implicit and explicit parts of a `TransientIMEXFEOperator` must be
  defined on the same test and trial spaces and have the same assembler.
  """
  @assert (get_test(im_tfeop) == get_test(ex_tfeop)) msg
  @assert (get_trial(im_tfeop) == get_trial(ex_tfeop)) msg
  @assert (get_assembler(im_tfeop) == get_assembler(ex_tfeop)) msg

  im_order, ex_order = get_order(im_tfeop), get_order(ex_tfeop)
  check_imex_compatibility(im_order, ex_order)
end

function IMEXODEOperatorType(
  im_tfeop::TransientFEOperator, ex_tfeop::TransientFEOperator
)
  T_im, T_ex = ODEOperatorType(im_tfeop), ODEOperatorType(ex_tfeop)
  IMEXODEOperatorType(T_im, T_ex)
end

##################################
# GenericTransientIMEXFEOperator #
##################################
"""
    struct GenericTransientIMEXFEOperator <: TransientIMEXFEOperator end
"""
struct GenericTransientIMEXFEOperator{T<:ODEOperatorType} <: TransientIMEXFEOperator{T}
  im_tfeop::TransientFEOperator
  ex_tfeop::TransientFEOperator

  function GenericTransientIMEXFEOperator(
    im_tfeop::TransientFEOperator,
    ex_tfeop::TransientFEOperator
  )
    check_imex_compatibility(im_tfeop, ex_tfeop)
    T = IMEXODEOperatorType(im_tfeop, ex_tfeop)
    new{T}(im_tfeop, ex_tfeop)
  end
end

# Default constructor
function TransientIMEXFEOperator(
  im_tfeop::TransientFEOperator,
  ex_tfeop::TransientFEOperator
)
  GenericTransientIMEXFEOperator(im_tfeop, ex_tfeop)
end

# TransientIMEXFEOperator interface
function get_imex_operators(tfeop::GenericTransientIMEXFEOperator)
  (tfeop.im_tfeop, tfeop.ex_tfeop)
end

########
# Test #
########
"""
    test_tfe_operator(
      tfeop::TransientFEOperator,
      t::Real, uh::TransientCellField
    ) -> Bool

Test the interface of `TransientFEOperator` specializations.
"""
function test_tfe_operator(
  tfeop::TransientFEOperator,
  t::Real, uh::TransientCellField
)
  U = get_trial(tfeop)
  Ut = U(t)
  @test Ut isa FESpace

  V = get_test(tfeop)
  @test V isa FESpace

  odeop = get_algebraic_operator(tfeop)
  @test odeop isa ODEOperator

  us = (get_free_dof_values(uh.cellfield),)
  for derivative in uh.derivatives
    us = (us..., get_free_dof_values(derivative))
  end

  test_ode_operator(odeop, t, us)

  true
end
