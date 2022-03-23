"""
It represents a FE function at a set of time steps. It is a wrapper of a ODE
solution for free values combined with data for Dirichlet values. Thus, it is a
lazy iterator that computes the solution at each time step when accessing the
solution.
"""
struct TransientFESolution
  odesol::ODESolution
  trial
end


function TransientFESolution(
  solver::ODESolver, op::TransientFEOperator, uh0, t0::Real, tF::Real)

  ode_op = get_algebraic_operator(op)
  u0 = get_free_dof_values(uh0)
  ode_sol = solve(solver,ode_op,u0,t0,tF)
  trial = get_trial(op)

  TransientFESolution(ode_sol, trial)
end

function TransientFESolution(
  solver::ODESolver,
  op::TransientFEOperator,
  xh0::Tuple{Vararg{Any}},
  t0::Real,
  tF::Real)

  ode_op = get_algebraic_operator(op)
  x0 = ()
  for xhi in xh0
    x0 = (x0...,get_free_dof_values(xhi))
  end
  ode_sol = solve(solver,ode_op,x0,t0,tF)
  trial = get_trial(op)

  TransientFESolution(ode_sol, trial)
end

# Solve functions

function solve(
  solver::ODESolver,op::TransientFEOperator,u0,t0::Real,tf::Real)
  TransientFESolution(solver,op,u0,t0,tf)
end

function solve(
  solver::ODESolver,op::TransientFEOperator,u0,v0,a0,t0::Real,tf::Real)
  TransientFESolution(solver,op,u0,v0,a0,t0,tf)
end

function test_transient_fe_solver(solver::ODESolver,op::TransientFEOperator,u0,t0,tf)
  solution = solve(solver,op,u0,t0,tf)
  test_transient_fe_solution(solution)
end

#@fverdugo this is a general implementation of iterate for TransientFESolution
# We could also implement another one for the very common case that the
# underlying ode_op is a ODEOpFromFEOp object

function Base.iterate(sol::TransientFESolution)

  odesolnext = Base.iterate(sol.odesol)

  if odesolnext === nothing
    return nothing
  end

  (uf, tf), odesolstate = odesolnext

  Uh = allocate_trial_space(sol.trial)
  Uh = evaluate!(Uh,sol.trial,tf)
  uh = FEFunction(Uh,uf)

  state = (Uh, odesolstate)

  (uh, tf), state
end

function Base.iterate(sol::TransientFESolution, state)

  Uh, odesolstate = state

  odesolnext = Base.iterate(sol.odesol,odesolstate)

  if odesolnext === nothing
    return nothing
  end

  (uf, tf), odesolstate = odesolnext

  Uh = evaluate!(Uh,sol.trial,tf)
  uh = FEFunction(Uh,uf)

  state = (Uh, odesolstate)

  (uh, tf), state

end

function test_transient_fe_solution(fesol::TransientFESolution)
  for (uhn,tn) in fesol
    @test isa(uhn,FEFunction)
    @test isa(tn,Real)
  end
  true
end
