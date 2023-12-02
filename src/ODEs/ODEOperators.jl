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
residual(t, us) = res(t, us[0], ..., us[N])
```, where `N` is the order of the ODE operator and `us[k] = ∂t^k(u)` is the
k-th-order time derivative of u.

# Mandatory
- [`get_order(odeop)`](@ref)
- [`allocate_residual(odeop, t, us, odeopcache)`](@ref)
- [`residual!(r, odeop, t, us, odeopcache)`](@ref)
- [`allocate_jacobian(odeop, t, us, odeopcache)`](@ref)
- [`jacobian!(J, odeop, t, us, k, γ, odeopcache)`](@ref)

# Optional
- [`allocate_odeopcache(odeop, args...)`](@ref)
- [`update_odeopcache!(odeopcache, odeop, t, args...)`](@ref)
- [`residual(odeop, t, us, odeopcache)`](@ref)
- [`jacobian(odeop, t, us, k, γ, odeopcache)`](@ref)
- [`jacobians!(J, odeop, t, us, γs, odeopcache)`](@ref)
- [`is_jacobian_constant(odeop, k)`](@ref)
- [`is_forcing_constant(odeop)`](@ref)
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
    allocate_odeopcache(odeop::ODEOperator, args...) -> CacheType

Allocate the cache required by the `ODEOperator`.
"""
function allocate_odeopcache(odeop::ODEOperator, args...)
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
    is_jacobian_constant(odeop::ODEOperator, k::Integer) -> Bool

Indicate whether the jacobian of the residual of the `ODEOperator` with respect
to the `k`-th-order time derivative is constant.
"""
function is_jacobian_constant(odeop::ODEOperator, k::Integer)
  false
end

"""
    is_forcing_constant(odeop::ODEOperator{<:AbstractMassLinearODE}) -> Bool

For an `ODEOperator` of type `AbstractMassLinearODE`, indicate whether the
forcing term is constant. For example with a `MassLinearODEOperator`,
```math
residual(t, u, v) = mass(t, ∂t(u), v) + res(t, u, v),```
this function indicates whether `res` is constant.
"""
function is_forcing_constant(odeop::ODEOperator{<:AbstractMassLinearODE})
  false
end

########
# Test #
########
"""
    test_ode_operator(
      odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      args...
    ) -> Bool

Test the interface of `ODEOperator` specializations.
"""
function test_ode_operator(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  args...
)
  odeopcache = allocate_odeopcache(odeop, args...)
  odeopcache = update_odeopcache!(odeopcache, odeop, t)

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
    @test is_jacobian_constant(odeop, k) isa Bool
  end

  if ODEOperatorType(odeop) <: AbstractMassLinearODE
    @test is_forcing_constant(odeop) isa Bool
  end

  true
end
