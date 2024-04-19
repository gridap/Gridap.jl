###############
# ODESolution #
###############
"""
    abstract type ODESolution <: GridapType end

Wrapper around an `ODEOperator` and `ODESolver` that represents the solution at
a set of time steps. It is an iterator that computes the solution at each time
step in a lazy fashion when accessing the solution.

# Mandatory
- [`iterate(odesltn)`](@ref)
- [`iterate(odesltn, state)`](@ref)
"""
abstract type ODESolution <: GridapType end

"""
    Base.iterate(odesltn::ODESolution) -> ((Real, AbstractVector), StateType)

Allocate the operators and cache and perform one time step of the `ODEOperator`
with the `ODESolver` attached to the `ODESolution`.
"""
function Base.iterate(odesltn::ODESolution)
  @abstractmethod
end

"""
    Base.iterate(odesltn::ODESolution) -> ((Real, AbstractVector), StateType)

Perform one time step of the `ODEOperator` with the `ODESolver` attached to the
`ODESolution`.
"""
function Base.iterate(odesltn::ODESolution, state)
  @abstractmethod
end

Base.IteratorSize(::Type{<:ODESolution}) = Base.SizeUnknown()

######################
# GenericODESolution #
######################
"""
    struct GenericODESolution <: ODESolution end

Generic wrapper for the evolution of an `ODEOperator` with an `ODESolver`.
"""
struct GenericODESolution <: ODESolution
  odeslvr::ODESolver
  odeop::ODEOperator
  t0::Real
  tF::Real
  us0::Tuple{Vararg{AbstractVector}}
end

function Base.iterate(odesltn::GenericODESolution)
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  t0, us0 = odesltn.t0, odesltn.us0

  # Allocate cache
  odecache = allocate_odecache(odeslvr, odeop, t0, us0)

  # Starting procedure
  state0, odecache = ode_start(
    odeslvr, odeop,
    t0, us0,
    odecache
  )

  # Marching procedure
  stateF = copy.(state0)
  tF, stateF, odecache = ode_march!(
    stateF,
    odeslvr, odeop,
    t0, state0,
    odecache
  )

  # Finishing procedure
  uF = copy(first(us0))
  uF, odecache = ode_finish!(
    uF,
    odeslvr, odeop,
    t0, tF, stateF,
    odecache
  )

  # Update iterator
  data = (tF, uF)
  state = (tF, stateF, state0, uF, odecache)
  (data, state)
end

function Base.iterate(odesltn::GenericODESolution, state)
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  t0, state0, stateF, uF, odecache = state

  if t0 >= odesltn.tF - ε
    return nothing
  end

  # Marching procedure
  tF, stateF, odecache = ode_march!(
    stateF,
    odeslvr, odeop,
    t0, state0,
    odecache
  )

  # Finishing procedure
  uF, odecache = ode_finish!(
    uF,
    odeslvr, odeop,
    t0, tF, stateF,
    odecache
  )

  # Update iterator
  data = (tF, uF)
  state = (tF, stateF, state0, uF, odecache)
  (data, state)
end

##############################
# Default behaviour of solve #
##############################
"""
    solve(
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, tF::Real, us0::Tuple{Vararg{AbstractVector}},
    ) -> ODESolution

Create an `ODESolution` wrapper around the `ODEOperator` and `ODESolver`,
starting with state `us0` at time `t0`, to be evolved until `tF`.
"""
function Algebra.solve(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, tF::Real, us0::Tuple{Vararg{AbstractVector}},
)
  GenericODESolution(odeslvr, odeop, t0, tF, us0)
end

function Algebra.solve(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, tF::Real, u0::AbstractVector,
)
  us0 = (u0,)
  solve(odeslvr, odeop, t0, tF, us0)
end

########
# Test #
########
"""
    test_ode_solution(odesltn::ODESolution) -> Bool

Test the interface of `ODESolution` specializations.
"""
function test_ode_solution(odesltn::ODESolution)
  for (t_n, us_n) in odesltn
    @test t_n isa Real
    @test us_n isa AbstractVector
  end
  true
end
