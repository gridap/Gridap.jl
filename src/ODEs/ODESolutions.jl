###############
# ODESolution #
###############
"""
    abstract type ODESolution <: GridapType end

Wrapper around an `ODEOperator` and an `ODESolver`. It is an iterator that
computes the solution at each time step in a lazy fashion when accessing the
solution.

# Mandatory
- [`Base.iterate(odesltn)`](@ref)
- [`Base.iterate(odesltn, state)`](@ref)
"""
abstract type ODESolution <: GridapType end

"""
    Base.iterate(odesltn::ODESolution) -> ((Tuple{Vararg{AbstractVector}}, Real), StateType)

Allocate a cache and perform one time step of the `ODEOperator` with the
`ODESolver` attached to the `ODESolution`.
"""
function Base.iterate(odesltn::ODESolution)
  @abstractmethod
end

"""
    Base.iterate(odesltn::ODESolution) -> ((Tuple{Vararg{AbstractVector}}, Real), StateType)

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
  us0::Tuple{Vararg{AbstractVector}}
  t0::Real
  tF::Real
end

function GenericODESolution(
  odeslvr::ODESolver, odeop::ODEOperator,
  us0::AbstractVector, t0::Real, tF::Real
)
  GenericODESolution(odeslvr, odeop, (us0,), t0, tF)
end

function Base.iterate(odesltn::GenericODESolution)
  us0 = ()
  usF = ()
  for k in eachindex(odesltn.us0)
    us0 = (us0..., copy(odesltn.us0[k]))
    usF = (usF..., copy(odesltn.us0[k]))
  end
  t0 = odesltn.t0

  # Solve step
  usF, tF, cache = solve_step!(usF, odesltn.odeslvr, odesltn.odeop, us0, t0)

  # Update state
  for k in eachindex(us0, usF)
    copy!(us0[k], usF[k])
  end
  state = (usF, us0, tF, cache)

  (first(usF), tF), state
end

function Base.iterate(odesltn::GenericODESolution, state)
  usF, us0, t0, cache = state

  if t0 >= odesltn.tF - ε
    return nothing
  end

  # Solve step
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  usF, tF, cache = solve_step!(usF, odeslvr, odeop, us0, t0, cache)

  # Update state
  for k in eachindex(us0, usF)
    copy!(us0[k], usF[k])
  end
  state = (usF, us0, tF, cache)

  (first(usF), tF), state
end

########
# Test #
########
"""
    test_ode_solution(odesltn::ODESolution) -> Bool

Test the interface of `ODESolution` specializations.
"""
function test_ode_solution(odesltn::ODESolution)
  for (us_n, t_n) in odesltn
    @test t_n isa Real
    @test us_n isa AbstractVector
  end
  true
end
