###################
# ODEOperatorType #
###################
"""
Trait for ODEOperator that indicates the type of ODE (nonlinear, linear with
respect to the highest time derivative, with constant coefficients).
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
Nonlinear ODE defined by a residual of the form
```math
residual(t, u, v) = res(t, u, ∂t(u), ..., ∂t^N(u), v)
```, where `∂t^k(u)` is the k-th time derivative of u and `N` is the order of
the ODE operator.
"""
abstract type ODEOperator{T<:ODEOperatorType} <: GridapType end

"""
ODE whose residual is linear with respect to the highest-order time derivative,
e.g.
```math
residual(t, u, v) = mass(t, ∂t^N(u), v) + res(t, u, ∂t(u), ..., ∂t^(N-1)(u), v)
```, where `N` is the order of the ODE operator, `mass` is linear in `∂t^N(u)`
and `res` does not depend on `∂t^N(u)`.
"""
const MassLinearODEOperator = ODEOperator{MassLinearODE}

"""
ODE whose residual is linear with respect to the highest-order time derivative,
and whose corresponding matrix is constant with respect to time and lower-order
time derivatives, e.g.
```math
residual(t, u, v) = mass(∂t^N(u), v) + res(t, u, ∂t(u), ..., ∂t^(N-1)(u), v)
```, where `N` is the order of the ODE operator, `mass` is linear in `∂t^N(u)`
and independent of time, and `res` does not depend on `∂t^N(u)`.
"""
const ConstantMassODEOperator = ODEOperator{ConstantMassODE}

"""
Return the `ODEOperatorType` of an ODE operator.
"""
ODEOperatorType(::ODEOperator{T}) where {T} = T

"""
Return the order of the ODE operator
"""
function get_order(::ODEOperator)
  @abstractmethod
end

"""
Allocate a residual vector for the ODE operator
"""
function allocate_residual(
  op::ODEOperator,
  t::Real, us::VecOrNTupleVec,
  cache
)
  @abstractmethod
end

"""
Return the residual vector of the ODE operator at a given point
(t, u, ∂t(u), ..., ∂t^n(u)), where `n` is the order of the ODE
"""
function residual!(
  r::AbstractVector, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  @abstractmethod
end

"""
Allocate a jacobian matrix for the ODE operator
"""
function allocate_jacobian(
  op::ODEOperator,
  t::Real, us::VecOrNTupleVec,
  cache
)
  @abstractmethod
end

"""
Add the jacobian with respect to the i-th time derivative, weighted by some
factor γ, to the matrix J.
"""
function jacobian!(
  J::AbstractMatrix, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  i::Integer, γ::Real,
  cache
)
  @abstractmethod
end

"""
Add a linear combination of all jacobians, weighted by some weights γs
"""
function jacobians!(
  J::AbstractMatrix, op::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  @abstractmethod
end

"""
Allocate the cache required by the `ODESolution` for the ODE operator
"""
function allocate_cache(op::ODEOperator, args...)
  @abstractmethod
end

"""
Update the cache of the `ODESolution` attached to the ODE operator
"""
function update_cache!(cache, op::ODEOperator, t::Real)
  @abstractmethod
end

########
# Test #
########
"""
Test the interface of `ODEOperator` specializations
"""
function test_ode_operator(
  op::ODEOperator,
  t::Real, u::AbstractVector, u̇::AbstractVector
)
  cache = allocate_cache(op)
  cache = update_cache!(cache, op, t)
  r = allocate_residual(op, t, (u, u̇), cache)
  residual!(r, op, t, (u, u̇), cache)
  J = allocate_jacobian(op, t, (u, u̇), cache)
  jacobian!(J, op, t, (u, u̇), 0, 1, cache)
  jacobian!(J, op, t, (u, u̇), 1, 1, cache)
  jacobians!(J, op, t, (u, u̇), (1, 1), cache)
  true
end
