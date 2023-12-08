###################
# ODEOperatorType #
###################
"""
    abstract type ODEOperatorType <: GridapType end

Trait that indicates the (linearity) type of an ODE operator.
"""
abstract type ODEOperatorType <: GridapType end
struct NonlinearODE <: ODEOperatorType end

abstract type AbstractQuasilinearODE <: ODEOperatorType end
struct QuasilinearODE <: AbstractQuasilinearODE end

abstract type AbstractSemilinearODE <: AbstractQuasilinearODE end
struct SemilinearODE <: AbstractSemilinearODE end

abstract type AbstractLinearODE <: AbstractSemilinearODE end
struct LinearODE <: AbstractLinearODE end

###############
# ODEOperator #
###############
"""
    abstract type ODEOperator <: GridapType end

General implicit, nonlinear ODE operator defined by a residual of the form
```math
residual(t, ∂t^0[u], ..., ∂t^N[u]) = 0,
```
where `N` is the order of the ODE operator and `∂t^k[u]` is the `k`-th-order
time derivative of `u`.

# Mandatory
- [`get_order(odeop)`](@ref)
- [`allocate_residual(odeop, t, us, odeopcache)`](@ref)
- [`residual!(r, odeop, t, us, odeopcache)`](@ref)
- [`allocate_jacobian(odeop, t, us, odeopcache)`](@ref)
- [`jacobian!(J, odeop, t, us, k, γ, odeopcache)`](@ref)

# Optional
- [`allocate_odeopcache(odeop, t, us, args...)`](@ref)
- [`update_odeopcache!(odeopcache, odeop, t, args...)`](@ref)
- [`residual(odeop, t, us, odeopcache)`](@ref)
- [`jacobian(odeop, t, us, k, γ, odeopcache)`](@ref)
- [`jacobians!(J, odeop, t, us, γs, odeopcache)`](@ref)
- [`is_form_constant(odeop, k)`](@ref)
"""
abstract type ODEOperator{C<:ODEOperatorType} <: GridapType end

"""
    QuasilinearODEOperator

ODE operator whose residual is linear with respect to the highest-order time
derivative, i.e.
```math
residual(t, ∂t^0[u], ..., ∂t^N[u]) = mass(t, ∂t^0[u], ..., ∂t^(N-1)[u]) ∂t^N[u]
                                   +  res(t, ∂t^0[u], ..., ∂t^(N-1)[u]),
```
where `N` is the order of the ODE operator, `∂t^k[u]` is the `k`-th-order time
derivative of `u`, and both `mass` and `res` have order `N-1`.

Alias for `ODEOperator{QuasilinearODE}`.
"""
const QuasilinearODEOperator = ODEOperator{QuasilinearODE}

"""
    SemilinearODEOperator

ODE operator whose residual is linear with respect to the highest-order time
derivative, and whose mass matrix only depend on time, i.e.
```math
residual(t, ∂t^0[u], ..., ∂t^N[u]) = mass(t) ∂t^N[u]
                                   +  res(t, ∂t^0[u], ..., ∂t^(N-1)[u]),
```
where `N` is the order of the ODE operator, `∂t^k[u]` is the `k`-th-order time
derivative of `u`, `mass` is independent of `u` and `res` has order `N-1`.

Alias for `ODEOperator{SemilinearODE}`.
"""
const SemilinearODEOperator = ODEOperator{SemilinearODE}

"""
    LinearODEOperator

ODE operator whose residual is linear with respect to all time derivatives, i.e.
```math
residual(t, ∂t^0[u], ..., ∂t^N[u]) = ∑_{0 ≤ k ≤ N} A_k(t) ∂t^k[u] + res(t),
```
where `N` is the order of the ODE operator, and `∂t^k[u]` is the `k`-th-order
time derivative of `u`.

Alias for `ODEOperator{LinearODE}`.
"""
const LinearODEOperator = ODEOperator{LinearODE}

"""
    ODEOperatorType(odeop::ODEOperator) -> ODEOperatorType

Return the `ODEOperatorType` of the `ODEOperator`.
"""
ODEOperatorType(odeop::ODEOperator) = ODEOperatorType(typeof(odeop))
ODEOperatorType(::Type{<:ODEOperator{C}}) where {C} = C

"""
    get_order(odeop::ODEOperator) -> Integer

Return the order of the `ODEOperator`.
"""
function Polynomials.get_order(odeop::ODEOperator)
  @abstractmethod
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
    residual(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache
    ) -> AbstractVector

Allocate a residual vector and evaluate it.
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
    residual!(
      r::AbstractVector, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache
    ) -> AbstractVector

Evaluate the residual vector of the `ODEOperator`.
"""
function Algebra.residual!(
  r::AbstractVector, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  @abstractmethod
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

"""
    jacobian(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      k::Integer, γ::Real,
      odeopcache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the `ODEOperator` and add the jacobian of the
residual of the `ODEOperator` with respect to the `k`-th-order time derivative,
weighted by some factor `γ`.
"""
function Algebra.jacobian(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, odeop, t, us, k, γ, odeopcache)
  J
end

"""
    jacobian!(
      J::AbstractMatrix, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      k::Integer, γ::Real,
      odeopcache
    ) -> AbstractMatrix

Add the jacobian of the residual of the `ODEOperator` with respect to the
`k`-th-order time derivative, weighted by some factor `γ`.
"""
function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  @abstractmethod
end

"""
    jacobians!(
      J::AbstractMatrix, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      γs::Tuple{Vararg{Real}},
      odeopcache
    ) -> AbstractMatrix

Add the jacobian of the residual of the `ODEOperator` with respect to all time
derivatives, weighted by some factors `γs`.
"""
function jacobians!(
  J::AbstractMatrix, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  for k in 0:get_order(odeop)
    γ = γs[k+1]
    if !iszero(γ)
      jacobian!(J, odeop, t, us, k, γ, odeopcache)
    end
  end
  J
end

"""
    is_form_constant(odeop::ODEOperator, k::Integer) -> Bool

Indicate whether the linear form of the `ODEOperator` corresponding to the
`k`-th-order time derivative of `u` is constant with respect to `t`.
"""
function is_form_constant(odeop::ODEOperator, k::Integer)
  false
end

###################
# IMEXODEOperator #
###################
"""
    abstract type IMEXODEOperator <: ODEOperator end

ODEOperator whose residual can be decomposed into
```math
residual(t, ∂t^0[u], ..., ∂t^N[u]) = implicit_residual(t, ∂t^0[u], ..., ∂t^N[u])
                                   + explicit_residual(t, ∂t^0[u], ..., ∂t^(N-1)[u]),
```
where
* The implicit operator defined by the implicit residual is considered stiff
and is meant to be solved with an implicit solver,
* The explicit operator defined by the explicit residual is considered non-stiff
and is meant to be solved with an explicit solver.

# Important
The explicit operator must have one order less than the implicit operator, so
that the mass term of the global operator is fully contained in the implicit
operator.

# Mandatory
- [`get_imex_operators(odeop)`](@ref)
"""
abstract type IMEXODEOperator{Cim,Cex} <: ODEOperator{Cim} end

"""
    get_imex_operators(odeop::IMEXODEOperator) -> (ODEOperator, ODEOperator)

Return the implicit and explicit `ODEOperator`s of the `IMEXODEOperator`.
"""
function get_imex_operators(odeop::IMEXODEOperator)
  @abstractmethod
end

# ODEOperator interface
function Polynomials.get_order(odeop::IMEXODEOperator)
  im_odeop, _ = get_imex_operators(odeop)
  get_order(im_odeop)
end

function allocate_odeopcache(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, args...
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache = allocate_odeopcache(im_odeop, t, im_us, args...)
  ex_odeopcache = allocate_odeopcache(ex_odeop, t, ex_us, args...)
  ex_res = allocate_residual(ex_odeop, t, ex_us, ex_odeopcache)
  (im_odeopcache, ex_odeopcache, ex_res)
end

function update_odeopcache!(
  odeopcache, odeop::IMEXODEOperator,
  t::Real, args...
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, ex_res = odeopcache
  update_odeopcache!(im_odeopcache, im_odeop, t, args...)
  update_odeopcache!(ex_odeopcache, ex_odeop, t, args...)
  (im_odeopcache, ex_odeopcache, ex_res)
end

function Algebra.allocate_residual(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, _ = odeopcache
  im_res = allocate_residual(im_odeop, t, im_us, im_odeopcache)
  ex_res = allocate_residual(ex_odeop, t, ex_us, ex_odeopcache)
  axpy!(1, ex_res, im_res)
  im_res
end

function Algebra.residual!(
  r::AbstractVector, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_res = r
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, ex_res = odeopcache
  residual!(im_res, im_odeop, t, im_us, im_odeopcache)
  residual!(ex_res, ex_odeop, t, ex_us, ex_odeopcache)
  axpy!(1, ex_res, im_res)
  r
end

function Algebra.allocate_jacobian(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, = odeopcache
  im_jac = allocate_jacobian(im_odeop, t, im_us, im_odeopcache)
  ex_jac = allocate_jacobian(ex_odeop, t, ex_us, ex_odeopcache)

  if !can_add_matrices(im_jac, ex_jac)
    throw("Cannot add matrices with different structures yet.")
  end
  axpy_sparse!(1, ex_jac, im_jac)
  im_jac
end

function Algebra.jacobian!(
  J::AbstractMatrix, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  im_us = us
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, = odeopcache
  jacobian!(J, im_odeop, t, im_us, k, γ, im_odeopcache)
  if k < get_order(odeop)
    ex_us = ntuple(i -> us[i], length(us) - 1)
    jacobian!(J, ex_odeop, t, ex_us, k, γ, ex_odeopcache)
  end
  J
end

function jacobians!(
  J::AbstractMatrix, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  im_us, ex_us = us, ntuple(i -> us[i], length(us) - 1)
  im_γs, ex_γs = γs, ntuple(i -> γs[i], length(γs) - 1)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, _ = odeopcache
  jacobians!(J, im_odeop, t, im_us, im_γs, im_odeopcache)
  jacobians!(J, ex_odeop, t, ex_us, ex_γs, ex_odeopcache)
  J
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

##########################
# GenericIMEXODEOperator #
##########################
"""
    struct GenericIMEXODEOperator <: IMEXODEOperator end

Generic `IMEXODEOperator`.
"""
struct GenericIMEXODEOperator{Cim,Cex} <: IMEXODEOperator{Cim,Cex}
  im_odeop::ODEOperator{Cim}
  ex_odeop::ODEOperator{Cex}

  function GenericIMEXODEOperator(im_odeop::ODEOperator, ex_odeop::ODEOperator)
    msg = """
    The explicit `ODEOperator` of an `IMEXODEOperator` must have one order less
    than the implicit `ODEOperator`.
    """
    @assert get_order(im_odeop) == get_order(ex_odeop) + 1 msg
    Cim = ODEOperatorType(im_odeop)
    Cex = ODEOperatorType(ex_odeop)
    new{Cim,Cex}(im_odeop, ex_odeop)
  end
end

function IMEXODEOperator(im_odeop::ODEOperator, ex_odeop::ODEOperator)
  GenericIMEXODEOperator(im_odeop, ex_odeop)
end

# IMEXODEOperator interface
function get_imex_operators(odeop::GenericIMEXODEOperator)
  (odeop.im_odeop, odeop.ex_odeop)
end

#########
# Utils #
#########
function can_add_matrices(A::AbstractMatrix, B::AbstractMatrix)
  if typeof(A) != typeof(B)
    return false
  end

  if A isa SparseMatrixCSC
    if rowvals(A) != rowvals(B)
      return false
    end

    if any(nzrange(A, j) != nzrange(B, j) for j in 1:size(A, 2))
      return false
    end
  end

  true
end

function axpy_sparse!(α::Real, A::AbstractMatrix, B::AbstractMatrix)
  axpy!(α, A, B)
end

function axpy_sparse!(α::Real, A::SparseMatrixCSC, B::SparseMatrixCSC)
  # TODO optimise the sum of sparse matrices

  # This is surprisingly better than axpy!(α, A, B)
  # @. B += α * A

  # This is way more efficient but only available when A and B have the same
  # structure (rowvals and nzrange)
  axpy!(α, nonzeros(A), nonzeros(B))
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
  odeopcache = allocate_odeopcache(odeop, t, us, args...)
  odeopcache = update_odeopcache!(odeopcache, odeop, t, args...)

  r = allocate_residual(odeop, t, us, odeopcache)
  @test r isa AbstractVector

  residual!(r, odeop, t, us, odeopcache)

  J = allocate_jacobian(odeop, t, us, odeopcache)
  @assert J isa AbstractMatrix

  order = get_order(odeop)
  for k in 0:order
    jacobian!(J, odeop, t, us, k, 1, odeopcache)
  end
  γs = ntuple(_ -> 1, order + 1)
  jacobians!(J, odeop, t, us, γs, odeopcache)

  for k in 0:order
    @test is_form_constant(odeop, k) isa Bool
  end

  true
end
