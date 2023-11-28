###############
# ODESolution #
###############
"""
    abstract type ODESolution <: GridapType end

Wrapper around an `ODEOperator` and an `ODESolver`. It is an iterator that
computes the solution at each time step in a lazy fashion when accessing the
solution.

# Mandatory
- [`Base.iterate(sol)`](@ref)
- [`Base.iterate(sol, state)`](@ref)
"""
abstract type ODESolution <: GridapType end

"""
    Base.iterate(sol::ODESolution) -> ((OneOrMoreVectors, Real), StateType)

Allocate a cache and perform one time step of the `ODEOperator` with the
`ODESolver` attached to the `ODESolution`.
"""
function Base.iterate(sol::ODESolution)
  @abstractmethod
end

"""
    Base.iterate(sol::ODESolution) -> ((OneOrMoreVectors, Real), StateType)

Perform one time step of the `ODEOperator` with the `ODESolver` attached to the
`ODESolution`.
"""
function Base.iterate(sol::ODESolution, state)
  @abstractmethod
end

Base.IteratorSize(::Type{<:ODESolution}) = Base.SizeUnknown()

######################
# GenericODESolution #
######################
"""
    struct GenericODESolution{T} <: ODESolution end

Generic wrapper for the evolution of an `ODEOperator` with an `ODESolver`.
"""
struct GenericODESolution{T} <: ODESolution
  solver::ODESolver
  op::ODEOperator
  us0::T
  t0::Real
  tF::Real
end

# Interface for AbstractVector-valued GenericODESolution
function Base.iterate(sol::GenericODESolution{<:AbstractVector})
  us0 = copy(sol.us0)
  usF = copy(sol.us0)
  t0 = sol.t0

  # Solve step
  usF, tF, cache = solve_step!(usF, sol.solver, sol.op, us0, t0)

  # Update state
  copy!(us0, usF)
  state = (usF, us0, tF, cache)

  (usF, tF), state
end

function Base.iterate(sol::GenericODESolution{<:AbstractVector}, state)
  usF, us0, t0, cache = state

  if t0 >= sol.tF - ε
    return nothing
  end

  # Solve step
  usF, tF, cache = solve_step!(usF, sol.solver, sol.op, us0, t0, cache)

  # Update state
  copy!(us0, usF)
  state = (usF, us0, tF, cache)

  (usF, tF), state
end

# Interface for Tuple{Vararg{AbstractVector}}-valued GenericODESolution
function Base.iterate(sol::GenericODESolution{<:Tuple{Vararg{AbstractVector}}})
  us0 = ()
  usF = ()
  for k in eachindex(sol.us0)
    us0 = (us0..., copy(sol.us0[k]))
    usF = (usF..., copy(sol.us0[k]))
  end
  t0 = sol.t0

  # Solve step
  usF, tF, cache = solve_step!(usF, sol.solver, sol.op, us0, t0)

  # Update state
  for k in eachindex(us0, usF)
    copy!(us0[k], usF[k])
  end
  state = (usF, us0, tF, cache)

  return (first(usF), tF), state
end

function Base.iterate(
  sol::GenericODESolution{<:Tuple{Vararg{AbstractVector}}},
  state
)
  usF, us0, t0, cache = state

  if t0 >= sol.tF - ε
    return nothing
  end

  # Solve step
  usF, tF, cache = solve_step!(usF, sol.solver, sol.op, us0, t0, cache)

  # Update state
  for k in eachindex(us0, usF)
    copy!(us0[k], usF[k])
  end
  state = (usF, us0, tF, cache)

  return (first(usF), tF), state
end

########
# Test #
########
"""
    test_ode_solution(sol::ODESolution) -> Bool

Test the interface of `ODESolution` specializations.
"""
function test_ode_solution(sol::ODESolution)
  for (us_n, t_n) in sol
    @test us_n isa AbstractVector
    @test t_n isa Real
  end
  true
end
