###################
# ODEOperatorType #
###################
"""
    abstract type ODEOperatorType <: GridapType end

Trait that indicates the (linearity) type of an ODE operator
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

General implicit, nonlinear ODE operator defined by a residual of the form
```math
residual(t, us) = res(t, us[1], ..., us[N])
```, where `N` is the order of the ODE operator and `us[k] = ∂t^k(u)` is the
k-th-order time derivative of u.

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

ODE operator whose residual is linear with respect to the highest-order time
derivative, e.g.
```math
residual(t, us) = mass(t, us[N]) + res(t, us[1], ..., us[N-1])
```, where `N` is the order of the ODE operator, `us[k] = ∂t^k(u)` is the
`k`-th-order time derivative of u, `mass` is linear in `∂t^N(u)` and `res` does
not depend on `∂t^N(u)`.

Alias for `ODEOperator{MassLinearODE}`.
"""
const MassLinearODEOperator = ODEOperator{MassLinearODE}

"""
    ConstantMassODEOperator

ODE operator whose residual is linear with respect to the highest-order time
derivative, and whose corresponding bilinear form is constant with respect to
time and lower-order time derivatives, e.g.
```math
residual(t, us) = mass(us[N]) + res(t, us[1], ..., us[N-1])
```, where `N` is the order of the ODE operator, `us[k] = ∂t^k(u)` is the
`k`-th-order time derivative of u, `mass` is linear in `∂t^N(u)` and independent
of time, and `res` does not depend on `∂t^N(u)`.

Alias for `ODEOperator{ConstantMassODE}`.
"""
const ConstantMassODEOperator = ODEOperator{ConstantMassODE}

"""
    ODEOperatorType(op::ODEOperator) -> ODEOperatorType

Return the `ODEOperatorType` of the `ODEOperator`
"""
ODEOperatorType(::ODEOperator{C}) where {C} = C

"""
    get_order(::ODEOperator) -> Integer

Return the order of the `ODEOperator`
"""
function Polynomials.get_order(::ODEOperator)
  @abstractmethod
end

"""
    allocate_cache(op::ODEOperator, args...) -> CacheType

Allocate the cache required by the `ODEOperator`
"""
function allocate_cache(op::ODEOperator, args...)
  @abstractmethod
end

"""
    update_cache!(cache, op::ODEOperator, t::Real, args...) -> CacheType

Update the cache of the `ODEOperator`
"""
function update_cache!(cache, op::ODEOperator, t::Real, args...)
  @abstractmethod
end

"""
    allocate_residual(
      op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      cache
    ) -> AbstractVector

Allocate a residual vector for the `ODEOperator`
"""
function Algebra.allocate_residual(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  @abstractmethod
end

"""
    residual(
      op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      cache
    ) -> AbstractVector

Allocate a residual vector and evaluate it
"""
function Algebra.residual(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  r = allocate_residual(op, t, us, cache)
  residual!(r, op, t, us, cache)
  r
end

"""
    residual!(
      r::AbstractVector, op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      cache
    ) -> AbstractVector

Evaluate the residual vector of the `ODEOperator`
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
      t::Real, us::Tuple{Vararg{AbstractVector}},
      cache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the `ODEOperator`
"""
function Algebra.allocate_jacobian(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  @abstractmethod
end

"""
    jacobian(
      op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      i::Integer, γ::Real,
      cache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the `ODEOperator` and add the jacobian of the
residual of the `ODEOperator` with respect to the `i`-th-order time derivative,
weighted by some factor `γ`
"""
function Algebra.jacobian(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  i::Integer, γ::Real,
  cache
)
  J = allocate_jacobian(op, t, us, cache)
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, op, t, us, i, γ, cache)
  J
end

"""
    jacobian!(
      J::AbstractMatrix, op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      i::Integer, γ::Real,
      cache
    ) -> AbstractMatrix

Add the jacobian of the residual of the `ODEOperator` with respect to the
`i`-th-order time derivative, weighted by some factor `γ`
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

Add the jacobian of the residual of the `ODEOperator` with respect to all time
derivatives, weighted by some factors `γs`
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
      t::Real, us::Tuple{Vararg{AbstractVector}}
    ) -> Bool

Test the interface of `ODEOperator` specializations
"""
function test_ode_operator(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}
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
