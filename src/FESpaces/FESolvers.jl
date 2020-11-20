
"""
"""
abstract type FESolver <: GridapType end

"""
    uh, cache = solve!(uh,solver,op)

This function changes the state of the input and can render it in a corrupted state.
It is recommended to rewrite the input `uh` with the output as illustrated to prevent any
issue.
"""
function solve!(uh::FEFunction,solver::FESolver,op::FEOperator)
  solve!(uh,solver,op,nothing)
end

"""
    uh, cache = solve!(uh,solver,op,cache)

This function changes the state of the input and can render it in a corrupted state.
It is recommended to rewrite the input `uh` with the output as illustrated to prevent any
issue. If `cache===nothing`, then it creates a new cache object.
"""
function solve!(uh::FEFunction,solver::FESolver,op::FEOperator,cache)
  @abstractmethod
end

"""
Solve that allocates, and sets initial guess to zero
and returns the solution
"""
function solve(nls::FESolver,op::FEOperator)
  U = get_trial(op)
  uh = zero(U)
  vh, cache = solve!(uh,nls,op)
  vh
end

function solve(op::AffineFEOperator)
  solver = LinearFESolver()
  solve(solver,op)
end

"""
"""
function solve(op::FEOperator)
  solver = NonlinearFESolver()
  solve(solver,op)
end

"""
"""
function test_fe_solver(
  nls::FESolver,
  op::FEOperator,
  x0::AbstractVector,
  x::AbstractVector,
  pred::Function=isapprox)

  trial = get_trial(op)

  u = FEFunction(trial,copy(x0))
  u, cache = solve!(u,nls,op)
  @test pred(get_free_values(u),x)

  u = FEFunction(trial,copy(x0))
  u,cache = solve!(u,nls,op,cache)
  @test pred(get_free_values(u),x)

  u = FEFunction(trial,copy(x0))
  u,cache = solve!(u,nls,op,cache)
  @test pred(get_free_values(u),x)

end

"""
The solver that solves a LinearFEOperator
"""
struct LinearFESolver <: FESolver
  ls::LinearSolver
end

"""
"""
function LinearFESolver()
  ls = LUSolver()
  LinearFESolver(ls)
end

function solve!(uh::FEFunction,solver::LinearFESolver,op::FEOperator, cache)
  @unreachable "Cannot solve a generic FEOperator with a LinearFESolver"
end

function solve!(u::FEFunction,solver::LinearFESolver,feop::AffineFEOperator,cache::Nothing)
  x = get_free_values(u)
  op = get_algebraic_operator(feop)
  cache = solve!(x,solver.ls,op)
  trial = get_trial(feop)
  u_new = FEFunction(trial,x)
  (u_new, cache)
end

function solve!(u::FEFunction,solver::LinearFESolver,feop::AffineFEOperator, cache)
  x = get_free_values(u)
  op = get_algebraic_operator(feop)
  cache = solve!(x,solver.ls,op,cache)
  trial = get_trial(feop)
  u_new = FEFunction(trial,x)
  (u_new,cache)
end

"""
A general NonlinearFESolver
"""
struct NonlinearFESolver <: FESolver
  nls::NonlinearSolver
end

"""
"""
function FESolver(nls::NonlinearSolver)
  NonlinearFESolver(nls)
end

function FESolver()
  NonlinearFESolver()
end

"""
"""
function NonlinearFESolver()
  nls = NLSolver(show_trace=false,method=:newton)
  NonlinearFESolver(nls)
end

function solve!(u::FEFunction,solver::NonlinearFESolver,feop::FEOperator,cache::Nothing)
  x = get_free_values(u)
  op = get_algebraic_operator(feop)
  cache = solve!(x,solver.nls,op)
  trial = get_trial(feop)
  u_new = FEFunction(trial,x)
  (u_new, cache)
end

function solve!(u::FEFunction,solver::NonlinearFESolver,feop::FEOperator,cache)
  x = get_free_values(u)
  op = get_algebraic_operator(feop)
  cache = solve!(x,solver.nls,op,cache)
  trial = get_trial(feop)
  u_new = FEFunction(trial,x)
  (u_new,cache)
end

