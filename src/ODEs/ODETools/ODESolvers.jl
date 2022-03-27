# Now, we need an abstract type representing a numerical discretization scheme
# for the ODE
"""
Represents a map that given (t_n,u_n) returns (t_n+1,u_n+1) and cache for the
corresponding `ODEOperator` and `NonlinearOperator`
"""
abstract type ODESolver <: GridapType end

function solve_step!(
  uF::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  solver::ODESolver,
  op::ODEOperator,
  u0::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  t0::Real,
  cache) # -> (uF,tF,cache)
  @abstractmethod
end

# Default API

function solve_step!(
  uF::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  solver::ODESolver,
  op::ODEOperator,
  u0::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  t0::Real) # -> (uF,tF,cache)
  solve_step!(uF,solver,op,u0,t0,nothing)
end

function solve(
  solver::ODESolver,
  op::ODEOperator,
  u0::T,
  t0::Real,
  tf::Real) where {T}
  GenericODESolution{T}(solver,op,u0,t0,tf)
end

# testers

function test_ode_solver(solver::ODESolver,op::ODEOperator,u0,t0,tf)
  solution = solve(solver,op,u0,t0,tf)
  test_ode_solution(solution)
end

# Specialization

include("ForwardEuler.jl")

include("ThetaMethod.jl")

include("AffineThetaMethod.jl")

include("RungeKutta.jl")

include("Newmark.jl")

include("AffineNewmark.jl")

include("ConstantNewmark.jl")

include("ConstantMatrixNewmark.jl")
