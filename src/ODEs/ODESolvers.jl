#######################
# DiscreteODEOperator #
#######################
"""
    abstract type DiscreteODEOperator <: NonlinearOperator end

Zero-th-order nonlinear operator (nonlinear system) corresponding to an
`ODEOperator` and an `ODESolver`.
"""
abstract type DiscreteODEOperator <: NonlinearOperator end

function Algebra.allocate_residual(
  disop::DiscreteODEOperator, x::AbstractVector
)
  zero(x)
end

#############################
# LinearDiscreteODEOperator #
#############################
"""
    abstract type LinearDiscreteODEOperator <: DiscreteODEOperator end

Zero-th-order linear operator (linear system) corresponding to an `ODEOperator`
and an `ODESolver`. This is mostly a copy of `AffineOperator` but enables the
interface with an `ODESolver`.

# Mandatory
- [`get_matrix(disop)`](@ref)
- [`get_vector(disop)`](@ref)
"""
abstract type LinearDiscreteODEOperator <: DiscreteODEOperator end

function Algebra.get_matrix(disop::LinearDiscreteODEOperator)
  @abstractmethod
end

function Algebra.get_vector(disop::LinearDiscreteODEOperator)
  @abstractmethod
end

function Algebra.residual!(
  r::AbstractVector,
  disop::LinearDiscreteODEOperator, x::AbstractVector
)
  mul!(r, get_matrix(disop), x)
  r .-= get_vector(disop)
  r
end

function Algebra.allocate_jacobian(
  disop::LinearDiscreteODEOperator, x::AbstractVector
)
  get_matrix(disop)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::LinearDiscreteODEOperator, x::AbstractVector
)
  copy_entries!(J, get_matrix(disop))
  J
end

#############
# ODESolver #
#############
"""
    abstract type ODESolver <: GridapType end

An `ODESolver` is a map that given `(t_n, us_n)` returns `(t_(n+1), us_(n+1))`
and the corresponding updated cache.

Here `us_n` is a vector of size `N-1`, where `N` is the order of the
`ODEOperator`, and `us_n[k] = ∂t^k(u)(t_n)` is the `k`-th-order time derivative
of `u` at `t_n`.

# Mandatory
- [`get_dt(odeslvr)`](@ref)
- [`allocate_disopcache(odeslvr, odeop, odeopcache, t, x)`](@ref)
- [`DiscreteODEOperator(odeslvr, odeop, odeopcache, disopcache, t, us, dt, args...)`](@ref)
- [`solve_step!(usF, odeslvr, odeop, us0, t0[, cache])`](@ref)

# Optional
- [`allocate_disslvrcache(odeslvr)`](@ref)
- [`solve(odeslvr, odeop, us0, t0, tF)`](@ref)
"""
abstract type ODESolver <: GridapType end

"""
    get_dt(odeslvr::ODESolver) -> Real

Return the time step of the `ODESolver`.
"""
function get_dt(odeslvr::ODESolver)
  @abstractmethod
end

"""
    allocate_disopcache(
      odeslvr::ODESolver,
      odeop::ODEOperator, odeopcache,
      t::Real, x::AbstractVector
    ) -> CacheType

Allocate the cache of the `DiscreteODEOperator` of the `ODESolver`.
"""
function allocate_disopcache(
  odeslvr::ODESolver,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  @abstractmethod
end

"""
    allocate_disslvrcache(odeslvr::ODESolver) -> CacheType

Allocate the cache of the solver for the `DiscreteODEOperator` corresponding to
the `ODESolver`.
"""
function allocate_disslvrcache(odeslvr::ODESolver)
  nothing
end

"""
    DiscreteODEOperator(
      odeslvr::ODESolver, odeop::ODEOperator,
      odeopcache, disopcache,
      t::Real, us::Tuple{Vararg{AbstractVector}}, dt::Real,
      args...
    ) -> DiscreteODEOperator

Return the `DiscreteODEOperator` corresponding to the `ODEOperator` and the
`ODESolver`.
"""
function DiscreteODEOperator(
  odeslvr::ODESolver, odeop::ODEOperator,
  odeopcache, disopcache,
  t::Real, us::Tuple{Vararg{AbstractVector}}, dt::Real,
  args...
)
  @abstractmethod
end

"""
    solve_step!(
      usF::Tuple{Vararg{AbstractVector}},
      odeslvr::ODESolver, odeop::ODEOperator,
      us0::Tuple{Vararg{AbstractVector}}, t0::Real
      [, cache]
    ) -> (Tuple{Vararg{AbstractVector}}, Real, CacheType)

Perform one time step of the `ODEOperator` with the `ODESolver` from `t0` with
initial state `us0`.
"""
function solve_step!(
  usF::Tuple{Vararg{AbstractVector}},
  odeslvr::ODESolver, odeop::ODEOperator,
  us0::Tuple{Vararg{AbstractVector}}, t0::Real,
  cache
)
  @abstractmethod
end

function solve_step!(
  usF::Tuple{Vararg{AbstractVector}},
  odeslvr::ODESolver, odeop::ODEOperator,
  us0::Tuple{Vararg{AbstractVector}}, t0::Real
)
  solve_step!(usF, odeslvr, odeop, us0, t0, nothing)
end

"""
    solve(
      odeslvr::ODESolver, odeop::ODEOperator,
      us0, t0::Real, tF::Real
    ) -> ODESolution

Create an `ODESolution` wrapper around the `ODEOperator` and `ODESolver`,
starting with state `us0` at time `t0`, to be evolved until `tF`.
"""
function Algebra.solve(
  odeslvr::ODESolver, odeop::ODEOperator,
  us0, t0::Real, tF::Real
)
  GenericODESolution(odeslvr, odeop, us0, t0, tF)
end

########
# Test #
########
"""
    test_ode_solver(
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, us0::Tuple{Vararg{AbstractVector}}, dt::Real, args...
    ) -> Bool

Test the interface of `ODESolver` specializations.
"""
function test_ode_solver(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}, dt::Real, args...
)
  @test get_dt(odeslvr) isa Real

  u0 = (us0 isa Tuple) ? first(us0) : us0
  odeopcache = allocate_odeopcache(odeop, t0, (us0..., u0))
  disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)

  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    t0, us0, dt,
    args...
  )
  @test disop isa DiscreteODEOperator

  usF = copy.(us0)
  (usF, tF, cache) = solve_step!(usF, odeslvr, odeop, us0, t0)
  (usF, tF, cache) = solve_step!(usF, odeslvr, odeop, us0, t0, cache)

  @test usF isa Tuple{Vararg{AbstractVector}}
  @test tF isa Real

  true
end

##################
# Import solvers #
##################
# First-order
include("ODESolvers/ForwardEuler.jl")

include("ODESolvers/ThetaMethod.jl")

include("ODESolvers/GeneralizedAlpha1.jl")

include("ODESolvers/Tableaus.jl")

include("ODESolvers/RungeKutta.jl")

include("ODESolvers/IMEXRungeKutta.jl")

# Second-order
include("ODESolvers/GeneralizedAlpha2.jl")

# TODO for now if a jacobian matrix is constant, it is not reassembled. This is
# nice, but we should also reduce the number of factorisations by storing them.
# This can always be done in the following scenarios
# * explicit schemes on `AbstractQuasilinearODE`s with constant mass
# * diagonally-implicit schemes on `LinearODE`s with constant mass and
# stiffness. Besides, it is often the case that different stages of a DIRK
# scheme have the same matrix. This happens when some diagonal coefficients have
# the same value. In the extreme case when all the diagonal values are the same,
# these methods are known as Singly-Diagonally-Implicit Runge-Kutta schemes
# (SDIRK). One special case happens when a diagonal coefficient is zero: the
# stage becomes explicit, even in the `AbstractQuasilinearODE` case. This
# strategy is already set up in the current implementation of DIRK.

# TODO another optimisation is the so-called FSAL property (First Same As Last)
# of some schemes, which can save one evaluation of the residual.
