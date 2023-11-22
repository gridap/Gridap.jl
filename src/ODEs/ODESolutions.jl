###############
# ODESolution #
###############
"""
    abstract type ODESolution <: GridapType end

Wrapper around an ODE operator and an ODE solution. It is an iterator that
computes the solution at each time step in a lazy fashion

- [`Base.iterate(sol)`](@ref)
- [`Base.iterate(sol, state)`](@ref)
"""
abstract type ODESolution <: GridapType end

"""
    Base.iterate(sol::ODESolution) -> ((OneOrMoreVectors, Real), StateType)

Allocate a cache and perform one time step of the ODE operator
"""
function Base.iterate(sol::ODESolution)
  @abstractmethod
end

"""
    Base.iterate(sol::ODESolution) -> ((OneOrMoreVectors, Real), StateType)

Perform one time step of the ODE operator
"""
function Base.iterate(sol::ODESolution, state)
  @abstractmethod
end

Base.IteratorSize(::Type{<:ODESolution}) = Base.SizeUnknown()

######################
# GenericODESolution #
######################
"""
    struct GenericODESolution{T} <: ODESolution

Generic wrapper for the evolution of an ODE operator with an ODE solver
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
  for i in eachindex(sol.us0)
    us0 = (us0..., copy(sol.us0[i]))
    usF = (usF..., copy(sol.us0[i]))
  end
  t0 = sol.t0

  # Solve step
  usF, tF, cache = solve_step!(usF, sol.solver, sol.op, us0, t0)

  # Update state
  for i in eachindex(us0, usF)
    copy!(us0[i], usF[i])
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
  for i in eachindex(us0, usF)
    copy!(us0[i], usF[i])
  end
  state = (usF, us0, tF, cache)

  return (first(usF), tF), state
end

########
# Test #
########
"""
    test_ode_solution(sol::ODESolution) -> Bool

Test the interface of `ODESolution` specializations
"""
function test_ode_solution(sol::ODESolution)
  for (us_n, t_n) in sol
    @test isa(us_n, AbstractVector)
    @test isa(t_n, Real)
  end
  true
end
