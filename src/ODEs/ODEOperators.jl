###################
# ODEOperatorType #
###################
"""
    abstract type ODEOperatorType <: GridapType end

Trait that indicates the linearity type of an ODE operator.
"""
abstract type ODEOperatorType <: GridapType end
struct NonlinearODE <: ODEOperatorType end

"""
    abstract type AbstractQuasilinearODE <: ODEOperatorType end

ODE operator whose residual is linear with respect to the highest-order time
derivative, i.e.
```
residual(t, ∂t^0[u], ..., ∂t^N[u]) = mass(t, ∂t^0[u], ..., ∂t^(N-1)[u]) ∂t^N[u]
                                   +  res(t, ∂t^0[u], ..., ∂t^(N-1)[u]),
```
where `N` is the order of the ODE operator, `∂t^k[u]` is the `k`-th-order time
derivative of `u`, and both `mass` and `res` have order `N-1`.
"""
abstract type AbstractQuasilinearODE <: ODEOperatorType end
struct QuasilinearODE <: AbstractQuasilinearODE end

"""
    abstract type AbstractSemilinearODE <: AbstractQuasilinearODE end

ODE operator whose residual is linear with respect to the highest-order time
derivative, and whose mass matrix only depend on time, i.e.
```
residual(t, ∂t^0[u], ..., ∂t^N[u]) = mass(t) ∂t^N[u]
                                   +  res(t, ∂t^0[u], ..., ∂t^(N-1)[u]),
```
where `N` is the order of the ODE operator, `∂t^k[u]` is the `k`-th-order time
derivative of `u`, `mass` is independent of `u` and `res` has order `N-1`.
"""
abstract type AbstractSemilinearODE <: AbstractQuasilinearODE end
struct SemilinearODE <: AbstractSemilinearODE end

"""
    abstract type AbstractLinearODE <: AbstractSemilinearODE end

ODE operator whose residual is linear with respect to all time derivatives, i.e.
```
residual(t, ∂t^0[u], ..., ∂t^N[u]) = ∑_{0 ≤ k ≤ N} A_k(t) ∂t^k[u] - res(t),
```
where `N` is the order of the ODE operator, and `∂t^k[u]` is the `k`-th-order
time derivative of `u`.
"""
abstract type AbstractLinearODE <: AbstractSemilinearODE end
struct LinearODE <: AbstractLinearODE end

################
# IMEX Helpers #
################
"""
    check_imex_compatibility(im_order::Integer, ex_order::Integer) -> Bool

Check whether two operators can make a valid IMEX operator decomposition. This
function should be called in the constructors of concrete IMEX operators.
"""
function check_imex_compatibility(im_order::Integer, ex_order::Integer)
  msg = """
  The explicit operator of an IMEX operator decomposition must have one order
  less than the implicit operator.
  """
  @assert (im_order == ex_order + 1) msg
end

"""
    IMEXODEOperatorType(
      T_im::Type{<:ODEOperatorType},
      T_ex::Type{<:ODEOperatorType}
    ) -> ODEOperatorType

Return the `ODEOperatorType` of the operator defined by an IMEX decomposition.
This function should be called in the constructors of concrete IMEX operators.
"""
function IMEXODEOperatorType(
  T_im::Type{<:ODEOperatorType},
  T_ex::Type{<:ODEOperatorType}
)
  T_im
end

function IMEXODEOperatorType(
  T_im::Type{<:AbstractLinearODE},
  T_ex::Type{<:ODEOperatorType}
)
  SemilinearODE
end

# We should theoretically dispatch on T_ex <: AbstractQuasilinearODE because
# in that case we can write the decomposition as
#   im_A_N(t) ∂t^N[u]
# + [im_A_(N-1)(t) + ex_mass(t, ∂t^0[u], ..., ∂t^(N-2)[u])] ∂t^(N-1)[u]
# + ∑_{0 ≤ k ≤ N-2} im_A_k(t) ∂t^k[u] + im_res(t) + ex_res(t, ∂t^0[u], ..., ∂t^(N-1)[u])
# so we can identify two linear forms corresponding to the two highest-order
# time derivatives, and then the rest of the residual. We decide to still
# define the global operator as semilinear for the following reasons:
# * For a first-order ODE, the explicit part has order zero, so the definitions
# of quasilinear, semilinear and linear coincide. This will default to the
# case below. This means there is only a special case when the residual has
# order two or higher.
# * We would need to have a new type when we can identify two linear forms
# corresponding to the two highest-order time derivatives. This would
# recursively force us to create order-dependent linearity types based on how
# many linear forms have been identified.
# * This distinction is not common in the litterature and indeed there does not
# seem to exist ODE solvers that take advantage of this kind of multi-form
# operator decomposition.

function IMEXODEOperatorType(
  T_im::Type{<:AbstractLinearODE},
  T_ex::Type{<:AbstractLinearODE}
)
  T_im
end

###############
# ODEOperator #
###############
"""
    abstract type ODEOperator <: GridapType end

General implicit, nonlinear ODE operator defined by a residual of the form
```
residual(t, ∂t^0[u], ..., ∂t^N[u]) = 0,
```
where `N` is the order of the ODE operator and `∂t^k[u]` is the `k`-th-order
time derivative of `u`.

# Mandatory
- [`get_order(odeop)`](@ref)
- [`get_forms(odeop)`](@ref)
- [`allocate_residual(odeop, t, us, odeopcache)`](@ref)
- [`residual!(r, odeop, t, us, odeopcache; add::Bool)`](@ref)
- [`allocate_jacobian(odeop, t, us, odeopcache)`](@ref)
- [`jacobian_add!(J, odeop, t, us, ws, odeopcache)`](@ref)

# Optional
- [`get_num_forms(odeop)`](@ref)
- [`is_form_constant(odeop, k)`](@ref)
- [`allocate_odeopcache(odeop, t, us)`](@ref)
- [`update_odeopcache!(odeopcache, odeop, t)`](@ref)
- [`residual(odeop, t, us, odeopcache)`](@ref)
- [`jacobian!(odeop, t, us, ws, odeopcache)`](@ref)
- [`jacobian(odeop, t, us, ws, odeopcache)`](@ref)
"""
abstract type ODEOperator{T<:ODEOperatorType} <: GridapType end

"""
    ODEOperatorType(odeop::ODEOperator) -> ODEOperatorType

Return the `ODEOperatorType` of the `ODEOperator`.
"""
ODEOperatorType(::ODEOperator{T}) where {T} = T
ODEOperatorType(::Type{<:ODEOperator{T}}) where {T} = T

"""
    get_order(odeop::ODEOperator) -> Integer

Return the order of the `ODEOperator`.
"""
function Polynomials.get_order(odeop::ODEOperator)
  @abstractmethod
end

"""
    get_num_forms(odeop::ODEOperator) -> Integer

Return the number of linear forms of the `ODEOperator`. See [`get_forms`](@ref)
"""
function get_num_forms(odeop::ODEOperator)
  0
end

function get_num_forms(odeop::ODEOperator{<:AbstractQuasilinearODE})
  1
end

function get_num_forms(odeop::ODEOperator{<:AbstractLinearODE})
  get_order(odeop) + 1
end

"""
    get_forms(odeop::ODEOperator) -> Tuple{Vararg{Function}}

Return the linear forms of the `ODEOperator`:
* For a general ODE operator, return an empty tuple,
* For a quasilinear ODE operator, return a tuple with the mass matrix,
* For a linear ODE operator, return all the linear forms.
"""
function get_forms(odeop::ODEOperator)
  ()
end

function get_forms(odeop::ODEOperator{<:AbstractQuasilinearODE})
  @abstractmethod
end

"""
    is_form_constant(odeop::ODEOperator, k::Integer) -> Bool

Indicate whether the linear form of the `ODEOperator` corresponding to the
`k`-th-order time derivative of `u` is constant with respect to `t`.
"""
function is_form_constant(odeop::ODEOperator, k::Integer)
  false
end

"""
    allocate_odeopcache(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}, args...
    ) -> CacheType

Allocate the cache required by the `ODEOperator`.
"""
function allocate_odeopcache(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, args...
)
  nothing
end

"""
    update_odeopcache!(odeopcache, odeop::ODEOperator, t::Real, args...) -> CacheType

Update the cache of the `ODEOperator`.
"""
function update_odeopcache!(odeopcache, odeop::ODEOperator, t::Real, args...)
  odeopcache
end

"""
    allocate_residual(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache
    ) -> AbstractVector

Allocate a residual vector for the `ODEOperator`.
"""
function Algebra.allocate_residual(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  @abstractmethod
end

"""
    residual!(
      r::AbstractVector, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache; add::Bool=false
    ) -> AbstractVector

Compute the residual of the `ODEOperator`. If `add` is true, this function adds
to `r` instead of erasing it.
"""
function Algebra.residual!(
  r::AbstractVector, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  @abstractmethod
end

"""
    residual(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache
    ) -> AbstractVector

Allocate a vector and evaluate the residual of the `ODEOperator`.
"""
function Algebra.residual(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  r = allocate_residual(odeop, t, us, odeopcache)
  residual!(r, odeop, t, us, odeopcache)
  r
end

"""
    allocate_jacobian(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the `ODEOperator`.
"""
function Algebra.allocate_jacobian(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  @abstractmethod
end

const jacobian_weights_order_msg = """
The weights are ordered by increasing order of time derivative, i.e. the first
weight corresponds to `∂residual / ∂u` and the last to
`∂residual / ∂(d^N u / dt^N)`.
"""

"""
    jacobian_add!(
      J::AbstractMatrix, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
      odeopcache
    ) -> AbstractMatrix

Add the jacobian of the residual of the `ODEOperator` with respect to all time
derivatives, weighted by some factors `ws`.

$(jacobian_weights_order_msg)
"""
function jacobian_add!(
  J::AbstractMatrix, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  @abstractmethod
end

"""
    jacobian!(
      J::AbstractMatrix, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
      odeopcache
    ) -> AbstractMatrix

Compute the jacobian of the residual of the `ODEOperator` with respect to all
time derivatives, weighted by some factors `ws`.

$(jacobian_weights_order_msg)
"""
function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  fillstored!(J, zero(eltype(J)))
  jacobian_add!(J, odeop, t, us, ws, odeopcache)
  J
end

"""
    jacobian(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
      odeopcache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the `ODEOperator` and compute the jacobian of
the residual of the `ODEOperator` with respect to all time derivatives,
weighted by some factors `ws`.

$(jacobian_weights_order_msg)
"""
function Algebra.jacobian(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  jacobian!(J, odeop, t, us, ws, odeopcache)
  J
end

###################
# IMEXODEOperator #
###################
"""
    abstract type IMEXODEOperator <: ODEOperator end

Implicit-Explicit decomposition of a residual defining an `ODEOperator`:
```
residual(t, ∂t^0[u], ..., ∂t^N[u]) = implicit_residual(t, ∂t^0[u], ..., ∂t^N[u])
                                   + explicit_residual(t, ∂t^0[u], ..., ∂t^(N-1)[u]),
```
where
* The implicit operator defined by the implicit residual is considered stiff and is meant to be solved implicitly,
* The explicit operator defined by the explicit residual is considered non-stiff and is meant to be solved explicitly.

# Important
The explicit operator must have one order less than the implicit operator, so
that the mass term of the global operator is fully contained in the implicit
operator.

# Mandatory
- [`get_imex_operators(odeop)`](@ref)
"""
abstract type IMEXODEOperator{T<:ODEOperatorType} <: ODEOperator{T} end

# IMEX Helpers
function check_imex_compatibility(im_odeop::ODEOperator, ex_odeop::ODEOperator)
  im_order, ex_order = get_order(im_odeop), get_order(ex_odeop)
  check_imex_compatibility(im_order, ex_order)
end

function IMEXODEOperatorType(im_odeop::ODEOperator, ex_odeop::ODEOperator)
  T_im, T_ex = ODEOperatorType(im_odeop), ODEOperatorType(ex_odeop)
  IMEXODEOperatorType(T_im, T_ex)
end

# IMEXODEOperator interface
"""
    get_imex_operators(odeop::IMEXODEOperator) -> (ODEOperator, ODEOperator)

Return the implicit and explicit parts of the `IMEXODEOperator`.
"""
function get_imex_operators(odeop::IMEXODEOperator)
  @abstractmethod
end

# ODEOperator interface
function Polynomials.get_order(odeop::IMEXODEOperator)
  im_odeop, _ = get_imex_operators(odeop)
  get_order(im_odeop)
end

function get_forms(odeop::IMEXODEOperator{<:AbstractQuasilinearODE})
  im_odeop, _ = get_imex_operators(odeop)
  im_forms = get_forms(im_odeop)
  (last(im_forms),)
end

function get_forms(odeop::IMEXODEOperator{<:AbstractLinearODE})
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_forms, ex_forms = get_forms(im_odeop), get_forms(ex_odeop)
  forms = ()
  for (im_form, ex_form) in zip(im_forms, ex_forms)
    form = (t, u) -> im_form(t, u) + ex_form(t, u)
    forms = (forms..., form)
  end
  (forms..., last(im_forms))
end

function is_form_constant(odeop::IMEXODEOperator, k::Integer)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_const = is_form_constant(im_odeop, k)
  ex_const = true
  if k < get_order(odeop)
    ex_const = is_form_constant(ex_odeop, k)
  end
  im_const && ex_const
end

function allocate_odeopcache(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, args...
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache = allocate_odeopcache(im_odeop, t, im_us, args...)
  ex_odeopcache = allocate_odeopcache(ex_odeop, t, ex_us, args...)
  (im_odeopcache, ex_odeopcache)
end

function update_odeopcache!(
  odeopcache, odeop::IMEXODEOperator,
  t::Real, args...
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache = odeopcache
  update_odeopcache!(im_odeopcache, im_odeop, t, args...)
  update_odeopcache!(ex_odeopcache, ex_odeop, t, args...)
  (im_odeopcache, ex_odeopcache)
end

function Algebra.allocate_residual(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache = odeopcache
  im_res = allocate_residual(im_odeop, t, im_us, im_odeopcache)
  ex_res = allocate_residual(ex_odeop, t, ex_us, ex_odeopcache)
  axpy!(1, ex_res, im_res)
  im_res
end

function Algebra.residual!(
  r::AbstractVector, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache = odeopcache
  residual!(r, im_odeop, t, im_us, im_odeopcache; add)
  residual!(r, ex_odeop, t, ex_us, ex_odeopcache; add=true)
  r
end

function Algebra.allocate_jacobian(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache = odeopcache

  # TODO Ideally, we want to allocate the jacobian matrix of both parts and sum them into
  # a new sparse matrix that has the sparsity structure of the sum. This is not fully
  # implemented for now.
  # * When both parts come from a TransientFEOperator, we replicate the code of
  # `allocate_jacobian` and simply merge the DomainContribution of both parts into a
  # single DomainContribution.
  # * Otherwise, for now, we allocate the two jacobians separately and add them. This will
  # break if they do not have the same sparsity structure.
  if im_odeop isa ODEOpFromTFEOp && ex_odeop isa ODEOpFromTFEOp
    # Common
    Ut = evaluate(get_trial(im_odeop.tfeop), nothing)
    du = get_trial_fe_basis(Ut)
    V = get_test(im_odeop.tfeop)
    v = get_fe_basis(V)
    assembler = get_assembler(im_odeop.tfeop)
    dc = DomainContribution()

    # Implicit part
    uh = _make_uh_from_us(im_odeop, us, im_odeopcache.Us)
    jacs = get_jacs(im_odeop.tfeop)
    for k in 0:get_order(im_odeop.tfeop)
      jac = jacs[k+1]
      dc = dc + jac(t, uh, du, v)
    end

    # Explicit part
    uh = _make_uh_from_us(ex_odeop, us, ex_odeopcache.Us)
    jacs = get_jacs(ex_odeop.tfeop)
    for k in 0:get_order(ex_odeop.tfeop)
      jac = jacs[k+1]
      dc = dc + jac(t, uh, du, v)
    end

    matdata = collect_cell_matrix(Ut, V, dc)
    allocate_matrix(assembler, matdata)
  else
    im_jac = allocate_jacobian(im_odeop, t, im_us, im_odeopcache)
    ex_jac = allocate_jacobian(ex_odeop, t, ex_us, ex_odeopcache)
    try
      axpy_entries!(1, ex_jac, im_jac)
    catch
      msg = """
      You are trying to define an IMEX operator where the jacobian of the implicit and
      explicit parts do not share the same sparsity structure. For now, this is only
      implemented when the implicit and explicit operators are `TransientFEOperator`.
      """
      @error msg
    end
    im_jac
  end
end

function jacobian_add!(
  J::AbstractMatrix, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_ws, ex_ws = ws, ntuple(i -> ws[i], length(ws) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache = odeopcache
  jacobian_add!(J, im_odeop, t, im_us, im_ws, im_odeopcache)
  jacobian_add!(J, ex_odeop, t, ex_us, ex_ws, ex_odeopcache)
  J
end

##########################
# GenericIMEXODEOperator #
##########################
"""
    struct GenericIMEXODEOperator <: IMEXODEOperator end

Generic `IMEXODEOperator`.
"""
struct GenericIMEXODEOperator{T} <: IMEXODEOperator{T}
  im_odeop::ODEOperator
  ex_odeop::ODEOperator

  function GenericIMEXODEOperator(im_odeop::ODEOperator, ex_odeop::ODEOperator)
    check_imex_compatibility(im_odeop, ex_odeop)
    T = IMEXODEOperatorType(im_odeop, ex_odeop)
    new{T}(im_odeop, ex_odeop)
  end
end

# Default constructor
function IMEXODEOperator(im_odeop::ODEOperator, ex_odeop::ODEOperator)
  GenericIMEXODEOperator(im_odeop, ex_odeop)
end

# IMEXODEOperator interface
function get_imex_operators(odeop::GenericIMEXODEOperator)
  (odeop.im_odeop, odeop.ex_odeop)
end

########
# Test #
########
"""
    test_ode_operator(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}, args...
    ) -> Bool

Test the interface of `ODEOperator` specializations.
"""
function test_ode_operator(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, args...
)
  num_forms = get_num_forms(odeop)
  for k in 0:num_forms-1
    @test is_form_constant(odeop, k) isa Bool
  end

  odeopcache = allocate_odeopcache(odeop, t, us, args...)
  odeopcache = update_odeopcache!(odeopcache, odeop, t, args...)

  r = allocate_residual(odeop, t, us, odeopcache)
  @test r isa AbstractVector

  residual!(r, odeop, t, us, odeopcache)

  J = allocate_jacobian(odeop, t, us, odeopcache)
  @assert J isa AbstractMatrix

  ws = ntuple(_ -> 1, get_order(odeop) + 1)
  jacobian!(J, odeop, t, us, ws, odeopcache)

  true
end
