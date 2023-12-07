###############
# ODESolution #
###############
"""
    abstract type ODESolution <: GridapType end

Wrapper around an `ODEOperator` and `ODESolver` that represents the solution at
a set of time steps. It is an iterator that computes the solution at each time
step in a lazy fashion when accessing the solution.

# Mandatory
- [`Base.iterate(odesltn)`](@ref)
- [`Base.iterate(odesltn, state)`](@ref)
"""
abstract type ODESolution <: GridapType end

"""
    Base.iterate(odesltn::ODESolution) -> ((Real, Tuple{Vararg{AbstractVector}}), StateType)

Allocate a cache and perform one time step of the `ODEOperator` with the
`ODESolver` attached to the `ODESolution`.
"""
function Base.iterate(odesltn::ODESolution)
  @abstractmethod
end

"""
    Base.iterate(odesltn::ODESolution) -> ((Real, Tuple{Vararg{AbstractVector}}), StateType)

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

function GenericODESolution(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, tF::Real, us0::AbstractVector,
)
  GenericODESolution(odeslvr, odeop, t0, tF, (us0,))
end

function Base.iterate(odesltn::GenericODESolution)
  t0 = odesltn.t0
  us0 = ()
  usF = ()
  for k in eachindex(odesltn.us0)
    ui0 = odesltn.us0[k]
    us0 = (us0..., copy(ui0))
    usF = (usF..., copy(ui0))
  end

  # Solve step
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  tF, usF, cache = solve_step!(usF, odeslvr, odeop, t0, us0)

  # Update state
  state = (tF, usF, us0, cache)

  (tF, first(usF)), state
end

function Base.iterate(odesltn::GenericODESolution, state)
  t0, us0, usF, cache = state

  if t0 >= odesltn.tF - ε
    return nothing
  end

  # Solve step
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  tF, usF, cache = solve_step!(usF, odeslvr, odeop, t0, us0, cache)

  # Update state
  state = (tF, usF, us0, cache)

  (tF, first(usF)), state
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
