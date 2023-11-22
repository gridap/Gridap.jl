###################
# ODEOperatorType #
###################
"""
    abstract type ODEOperatorType <: GridapType end

Trait that indicates the (linearity) type of an ODEOperator
"""
abstract type ODEOperatorType <: GridapType end
struct NonlinearODE <: ODEOperatorType end

abstract type AbstractMassLinearODE <: ODEOperatorType end
struct MassLinearODE <: AbstractMassLinearODE end
struct ConstantMassODE <: AbstractMassLinearODE end

###############
# ODEOperator #
###############
"""
    abstract type ODEOperator <: GridapType end

General implicit, nonlinear ODE defined by a residual of the form
```math
residual(t, u, v) = res(t, u, ∂t(u), ..., ∂t^N(u), v)
```, where `N` is the order of the ODE operator `∂t^k(u)` is the k-th time
derivative of u.

- [`get_order(op)`](@ref)
- [`allocate_cache(op, args...)`](@ref)
- [`update_cache!(cache, op, t, args...)`](@ref)
- [`allocate_residual(op, t, us, cache)`](@ref)
- [`residual!(r, op, t, us, cache)`](@ref)
- [`allocate_jacobian(op, t, us, cache)`](@ref)
- [`jacobian!(J, op, t, us, i, γ, cache)`](@ref)
- [`jacobians!(J, op, t, us, γs, cache)`](@ref)
"""
abstract type ODEOperator{C<:ODEOperatorType} <: GridapType end

"""
    MassLinearODEOperator

ODE whose residual is linear with respect to the highest-order time derivative,
e.g.
```math
residual(t, u, v) = mass(t, ∂t^N(u), v) + res(t, u, ∂t(u), ..., ∂t^(N-1)(u), v)
```, where `N` is the order of the ODE operator, `mass` is linear in `∂t^N(u)`
and `res` does not depend on `∂t^N(u)`.

Alias for `ODEOperator{MassLinearODE}`.
"""
const MassLinearODEOperator = ODEOperator{MassLinearODE}

"""
    ConstantMassODEOperator

ODE whose residual is linear with respect to the highest-order time derivative,
and whose corresponding bilinear form is constant with respect to time and
lower-order time derivatives, e.g.
```math
residual(t, u, v) = mass(∂t^N(u), v) + res(t, u, ∂t(u), ..., ∂t^(N-1)(u), v)
```, where `N` is the order of the ODE operator, `mass` is linear in `∂t^N(u)`
and independent of time, and `res` does not depend on `∂t^N(u)`.

Alias for `ODEOperator{ConstantMassODE}`.
"""
const ConstantMassODEOperator = ODEOperator{ConstantMassODE}

"""
    ODEOperatorType(op::ODEOperator) -> ODEOperatorType

Return the `ODEOperatorType` of the ODE operator
"""
ODEOperatorType(::ODEOperator{C}) where {C} = C

"""
    get_order(::ODEOperator) -> Integer

Return the order of the ODE operator
"""
function Polynomials.get_order(::ODEOperator)
  @abstractmethod
end

"""
    allocate_cache(op::ODEOperator, args...) -> CacheType

Allocate the cache required by the ODE operator
"""
function allocate_cache(op::ODEOperator, args...)
  @abstractmethod
end

"""
    update_cache!(cache, op::ODEOperator, t::Real, args...) -> CacheType

Update the cache of the ODE operator
"""
function update_cache!(cache, op::ODEOperator, t::Real, args...)
  @abstractmethod
end

"""
    allocate_residual(
      op::ODEOperator,
      t::Real, us::OneOrMoreVectors,
      cache
    ) -> AbstractVector

Allocate a residual vector for the ODE operator
"""
function Algebra.allocate_residual(
  op::ODEOperator,
  t::Real, us::OneOrMoreVectors,
  cache
)
  @abstractmethod
end

"""
    residual!(
      r::AbstractVector, op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      cache
    ) -> AbstractVector

Return the residual vector of the ODE operator at a given point
(t, u, ∂t(u), ..., ∂t^N(u)), where `N` is the order of the ODE
"""
function Algebra.residual!(
  r::AbstractVector, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  @abstractmethod
end

"""
    allocate_jacobian(
      op::ODEOperator,
      t::Real, us::OneOrMoreVectors,
      cache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the ODE operator
"""
function Algebra.allocate_jacobian(
  op::ODEOperator,
  t::Real, us::OneOrMoreVectors,
  cache
)
  @abstractmethod
end

"""
    jacobian!(
      J::AbstractMatrix, op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      i::Integer, γ::Real,
      cache
    ) -> AbstractMatrix

Add the jacobian with respect to the `i`-th time derivative, weighted by some
factor `γ`, to the matrix `J`
"""
function Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  i::Integer, γ::Real,
  cache
)
  @abstractmethod
end

"""
    jacobians!(
      J::AbstractMatrix, op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      γs::Tuple{Vararg{Real}},
      cache
    ) -> AbstractMatrix

Add the jacobian with respect to all time derivatives, weighted by some factors
`γs`, to the matrix `J`
"""
function jacobians!(
  J::AbstractMatrix, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  @abstractmethod
end

########
# Test #
########
"""
    test_ode_operator(
      op::ODEOperator,
      t::Real, us::OneOrMoreVectors
    ) -> Bool

Test the interface of `ODEOperator` specializations
"""
function test_ode_operator(
  op::ODEOperator,
  t::Real, us::OneOrMoreVectors
)
  cache = allocate_cache(op)
  cache = update_cache!(cache, op, t)

  r = allocate_residual(op, t, us, cache)
  residual!(r, op, t, us, cache)

  J = allocate_jacobian(op, t, us, cache)
  jacobian!(J, op, t, us, 0, 1, cache)
  jacobian!(J, op, t, us, 1, 1, cache)
  jacobians!(J, op, t, us, (1, 1), cache)

  true
end
