###################
# ODEOperatorType #
###################
"""
    abstract type ODEOperatorType <: GridapType end

Trait that indicates the (linearity) type of an ODE operator.
"""
abstract type ODEOperatorType <: GridapType end
struct NonlinearODE <: ODEOperatorType end

abstract type AbstractMassLinearODE <: ODEOperatorType end
struct MassLinearODE <: AbstractMassLinearODE end
struct LinearODE <: AbstractMassLinearODE end

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

# Mandatory
- [`get_order(op)`](@ref)
- [`allocate_residual(op, t, us, cache)`](@ref)
- [`residual!(r, op, t, us, cache)`](@ref)
- [`allocate_jacobian(op, t, us, cache)`](@ref)
- [`jacobian!(J, op, t, us, k, γ, cache)`](@ref)

# Optional
- [`allocate_cache(op, args...)`](@ref)
- [`update_cache!(cache, op, t, args...)`](@ref)
- [`residual(op, t, us, cache)`](@ref)
- [`jacobian(op, t, us, k, γ, cache)`](@ref)
- [`jacobians!(J, op, t, us, γs, cache)`](@ref)
- [`is_jacobian_constant(op, k)`](@ref)
- [`is_forcing_constant(op)`](@ref)
"""
abstract type ODEOperator{C<:ODEOperatorType} <: GridapType end

"""
    MassLinearODEOperator

ODE operator whose residual is linear with respect to the highest-order time
derivative, e.g.
```math
residual(t, us) = mass(t, us[1], ..., us[N-1]) us[N] + res(t, us[1], ..., us[N-1])
```, where `N` is the order of the ODE operator, `us[k] = ∂t^k(u)` is the
`k`-th-order time derivative of u, and `res` does not depend on `∂t^N(u)`.

Alias for `ODEOperator{MassLinearODE}`.
"""
const MassLinearODEOperator = ODEOperator{MassLinearODE}

"""
    LinearODEOperator

ODE operator whose residual is linear with respect to all time derivatives, e.g.
```math
residual(t, us) = A_N(t) us[N] + ... + A_1(t) us[1] + A_0(t) us[0] + res(t)
```, where `N` is the order of the ODE operator, `us[k] = ∂t^k(u)` is the
`k`-th-order time derivative of u.

Alias for `ODEOperator{LinearODE}`.
"""
const LinearODEOperator = ODEOperator{LinearODE}

"""
    ODEOperatorType(op::ODEOperator) -> ODEOperatorType

Return the `ODEOperatorType` of the `ODEOperator`.
"""
ODEOperatorType(op::ODEOperator) = ODEOperatorType(typeof(op))
ODEOperatorType(::Type{<:ODEOperator{C}}) where {C} = C

"""
    get_order(::ODEOperator) -> Integer

Return the order of the `ODEOperator`.
"""
function Polynomials.get_order(::ODEOperator)
  @abstractmethod
end

"""
    allocate_cache(op::ODEOperator, args...) -> CacheType

Allocate the cache required by the `ODEOperator`.
"""
function allocate_cache(op::ODEOperator, args...)
  nothing
end

"""
    update_cache!(cache, op::ODEOperator, t::Real, args...) -> CacheType

Update the cache of the `ODEOperator`.
"""
function update_cache!(cache, op::ODEOperator, t::Real, args...)
  cache
end

"""
    allocate_residual(
      op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      cache
    ) -> AbstractVector

Allocate a residual vector for the `ODEOperator`.
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

Allocate a residual vector and evaluate it.
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

Evaluate the residual vector of the `ODEOperator`.
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

Allocate a jacobian matrix for the `ODEOperator`.
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
      k::Integer, γ::Real,
      cache
    ) -> AbstractMatrix

Allocate a jacobian matrix for the `ODEOperator` and add the jacobian of the
residual of the `ODEOperator` with respect to the `k`-th-order time derivative,
weighted by some factor `γ`.
"""
function Algebra.jacobian(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  cache
)
  J = allocate_jacobian(op, t, us, cache)
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, op, t, us, k, γ, cache)
  J
end

"""
    jacobian!(
      J::AbstractMatrix, op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      k::Integer, γ::Real,
      cache
    ) -> AbstractMatrix

Add the jacobian of the residual of the `ODEOperator` with respect to the
`k`-th-order time derivative, weighted by some factor `γ`.
"""
function Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
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
derivatives, weighted by some factors `γs`.
"""
function jacobians!(
  J::AbstractMatrix, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  for k in 0:get_order(op)
    γ = γs[k+1]
    if !iszero(γ)
      jacobian!(J, op, t, us, k, γ, cache)
    end
  end
  J
end

"""
    is_jacobian_constant(op::ODEOperator, k::Integer) -> Bool

Indicate whether the jacobian of the residual of the `ODEOperator` with respect
to the `k`-th-order time derivative is constant.
"""
function is_jacobian_constant(op::ODEOperator, k::Integer)
  false
end

"""
    is_forcing_constant(op::ODEOperator{<:AbstractMassLinearODE}) -> Bool

For an `ODEOperator` of type `AbstractMassLinearODE`, indicate whether the
forcing term is constant. For example with a `MassLinearODEOperator`,
```math
residual(t, u, v) = mass(t, ∂t(u), v) + res(t, u, v),```
this function indicates whether `res` is constant.
"""
function is_forcing_constant(op::ODEOperator{<:AbstractMassLinearODE})
  false
end

########
# Test #
########
"""
    test_ode_operator(
      op::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}}
    ) -> Bool

Test the interface of `ODEOperator` specializations.
"""
function test_ode_operator(
  op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, args...
)
  cache = allocate_cache(op, args...)
  cache = update_cache!(cache, op, t)

  r = allocate_residual(op, t, us, cache)
  @test r isa AbstractVector

  residual!(r, op, t, us, cache)

  J = allocate_jacobian(op, t, us, cache)
  @assert J isa AbstractMatrix

  jacobian!(J, op, t, us, 0, 1, cache)
  jacobian!(J, op, t, us, 1, 1, cache)
  jacobians!(J, op, t, us, (1, 1), cache)

  for k in 0:get_order(op)
    @test is_jacobian_constant(op, k) isa Bool
  end

  if ODEOperatorType(op) <: AbstractMassLinearODE
    @test is_forcing_constant(op) isa Bool
  end

  true
end
