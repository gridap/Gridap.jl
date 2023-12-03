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
- [`residual!(r, odeop, t, us, odeopcache; include_mass)`](@ref)
- [`allocate_jacobian(odeop, t, us, odeopcache)`](@ref)
- [`jacobian!(J, odeop, t, us, k, γ, odeopcache)`](@ref)

# Optional
- [`allocate_odeopcache(odeop, t, us, args...)`](@ref)
- [`update_odeopcache!(odeopcache, odeop, t, args...)`](@ref)
- [`residual(odeop, t, us, odeopcache; include_mass)`](@ref)
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
      odeopcache; include_mass::Bool=true
    ) -> AbstractVector

Allocate a residual vector and evaluate it.
"""
function Algebra.residual(
  odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; include_mass::Bool=true
)
  r = allocate_residual(odeop, t, us, odeopcache)
  residual!(r, odeop, t, us, odeopcache; include_mass)
  r
end

"""
    residual!(
      r::AbstractVector, odeop::ODEOperator,
      t::Real, us::Tuple{Vararg{AbstractVector}},
      odeopcache; include_mass::Bool=true
    ) -> AbstractVector

Evaluate the residual vector of the `ODEOperator`.
"""
function Algebra.residual!(
  r::AbstractVector, odeop::ODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; include_mass::Bool=true
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
    is_forcing_constant(odeop::ODEOperator) -> Bool

For an `ODEOperator` of type `AbstractMassLinearODE`, indicate whether the
forcing term is constant. For example with a `MassLinearODEOperator`,
```math
residual(t, u, v) = mass(t, ∂t(u), v) + res(t, u, v),```
this function indicates whether `res` is constant.
"""
function is_forcing_constant(odeop::ODEOperator)
  false
end

###################
# IMEXODEOperator #
###################
"""
    abstract type IMEXODEOperator <: ODEOperator end

Pair of ODE operators: one operator is considered stiff and and meant to be
solved with an implicit solver, while the other is considered non stiff and
solved with an explicit solver.

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
  im_odeop, ex_odeop = get_imex_operators(odeop)
  get_order(im_odeop)
end

function allocate_odeopcache(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, args...
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache = allocate_odeopcache(im_odeop, t, us, args...)
  ex_cache = allocate_odeopcache(ex_odeop, t, us, args...)
  res_temp = allocate_residual(im_odeop, t, us, im_cache)
  (im_cache, ex_cache, res_temp)
end

function update_odeopcache!(
  odeopcache,
  odeop::IMEXODEOperator, t::Real, args...
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache, ex_cache, res_temp = odeopcache
  update_odeopcache!(im_cache, im_odeop, t, args...)
  update_odeopcache!(ex_cache, ex_odeop, t, args...)
  (im_cache, ex_cache, res_temp)
end

function Algebra.allocate_residual(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache, ex_cache, = odeopcache
  im_res = allocate_residual(im_odeop, t, us, im_cache)
  ex_res = allocate_residual(ex_odeop, t, us, ex_cache)
  axpy!(1, ex_res, im_res)
  im_res
end

function Algebra.residual!(
  r::AbstractVector, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; include_mass::Bool=true
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache, ex_cache, res_temp = odeopcache
  residual!(res_temp, im_odeop, t, us, im_cache; include_mass)
  residual!(r, ex_odeop, t, us, ex_cache; include_mass)
  axpy!(1, res_temp, r)
  r
end

function Algebra.allocate_jacobian(
  odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache, ex_cache, = odeopcache
  im_jac = allocate_jacobian(im_odeop, t, us, im_cache)
  ex_jac = allocate_jacobian(ex_odeop, t, us, ex_cache)

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
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache, ex_cache, = odeopcache
  jacobian!(J, im_odeop, t, us, k, γ, im_cache)
  jacobian!(J, ex_odeop, t, us, k, γ, ex_cache)
  J
end

function jacobians!(
  J::AbstractMatrix, odeop::IMEXODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_cache, ex_cache, = odeopcache
  jacobians!(J, im_odeop, t, us, γs, im_cache)
  jacobians!(J, ex_odeop, t, us, γs, ex_cache)
  J
end

function is_jacobian_constant(odeop::IMEXODEOperator, k::Integer)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_const = is_jacobian_constant(im_odeop, k)
  ex_const = is_jacobian_constant(ex_odeop, k)
  im_const && ex_const
end

function is_forcing_constant(odeop::IMEXODEOperator)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_const = is_forcing_constant(im_odeop)
  ex_const = is_forcing_constant(ex_odeop)
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
    An `IMEXODEOperator` can only be built from two `ODEOperator`s with same order.
    """
    @assert get_order(im_odeop) == get_order(ex_odeop) msg
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
  odeopcache = allocate_odeopcache(odeop, t, us, args...)
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
