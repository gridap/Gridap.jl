###############
# ODESolution #
###############
"""
Solution of an ODE during a given time interval. It is a lazy iterator that
computes the solution at each time step while it is iterated over
"""
abstract type ODESolution <: GridapType end

function iterate(sol::ODESolution)
  @abstractmethod
end

function iterate(sol::ODESolution, state)
  @abstractmethod
end

######################
# GenericODESolution #
######################
"""
Generic representation of the solution to an ODE problem
"""
struct GenericODESolution{T} <: ODESolution
  solver::ODESolver
  op::ODEOperator
  u0::T
  t0::Real
  tF::Real
end

# AbstractVector
function iterate(sol::GenericODESolution{<:AbstractVector})
  u0 = copy(sol.u0)
  uF = copy(sol.u0)
  t0 = sol.t0

  # Solve step
  uF, tF, cache = solve_step!(uF, sol.solver, sol.op, u0, t0)

  # Update state
  @. u0 = uF
  state = (uF, u0, tF, cache)

  (uF, tF), state
end

function iterate(sol::GenericODESolution{<:AbstractVector}, state)
  uF, u0, t0, cache = state

  if t0 >= sol.tF - ϵ
    return nothing
  end

  # Solve step
  uF, tF, cache = solve_step!(uF, sol.solver, sol.op, u0, t0, cache)

  # Update state
  @. u0 = uF
  state = (uF, u0, tF, cache)

  (uF, tF), state
end

# Tuple{Vararg{AbstractVector}}
function Base.iterate(sol::GenericODESolution{<:Tuple{Vararg{AbstractVector}}})
  u0 = ()
  uF = ()
  for i in eachindex(sol.u0)
    u0 = (u0..., copy(sol.u0[i]))
    uF = (uF..., copy(sol.u0[i]))
  end
  t0 = sol.t0

  # Solve step
  uF, tF, cache = solve_step!(uF, sol.solver, sol.op, u0, t0)

  # Update state
  for i in eachindex(u0, uF)
    @. u0[i] = uF[i]
  end
  state = (uF, u0, tF, cache)

  return (first(uF), tF), state
end

function Base.iterate(
  sol::GenericODESolution{<:Tuple{Vararg{AbstractVector}}},
  state
)
  uF, u0, t0, cache = state

  if t0 >= sol.tF - ϵ
    return nothing
  end

  # Solve step
  uF, tF, cache = solve_step!(uF, sol.solver, sol.op, u0, t0, cache)

  # Update state
  for i in eachindex(u0, uF)
    @. u0[i] = uF[i]
  end
  state = (uF, u0, tF, cache)

  return (first(uF), tF), state
end

########
# Test #
########
function test_ode_solution(sol::ODESolution)
  for (u_n, t_n) in sol
    @test isa(u_n, AbstractVector)
    @test isa(t_n, Real)
  end
  true
end
