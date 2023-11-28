#######################
# TransientFESolution #
#######################
"""
    struct TransientFESolution <: GridapType end

This represents an `FEFunction` at a set of time steps. It is a wrapper of an
`ODESolution` for free values combined with data for Dirichlet values. It is a
lazy iterator that computes the solution at each time step when accessing the
solution.
"""
struct TransientFESolution <: GridapType
  odesol::ODESolution
  trial
end

# Default constructors
function TransientFESolution(
  solver::ODESolver, fe_op::TransientFEOperator,
  uh0, t0::Real, tF::Real
)
  u0 = get_free_dof_values(uh0)
  ode_op = get_algebraic_operator(fe_op)
  ode_sol = solve(solver, ode_op, u0, t0, tF)
  trial = get_trial(fe_op)
  TransientFESolution(ode_sol, trial)
end

function TransientFESolution(
  solver::ODESolver, fe_op::TransientFEOperator,
  uh0::Tuple{Vararg{Any}}, t0::Real, tF::Real
)
  u0 = ()
  for uhi in uh0
    u0 = (u0..., get_free_dof_values(uhi))
  end
  ode_op = get_algebraic_operator(fe_op)
  ode_sol = solve(solver, ode_op, u0, t0, tF)
  trial = get_trial(fe_op)
  TransientFESolution(ode_sol, trial)
end

# ODESolver interface
function Algebra.solve(
  solver::ODESolver, op::TransientFEOperator,
  u0, t0::Real, tF::Real
)
  TransientFESolution(solver, op, u0, t0, tF)
end

# @fverdugo this is a general implementation of iterate for TransientFESolution
# We could also implement another one for the very common case that the
# underlying ode_op is a ODEOpFromFEOp object
function Base.iterate(sol::TransientFESolution)
  odesolnext = iterate(sol.odesol)
  if isnothing(odesolnext)
    return nothing
  end

  (uF, tF), odesolstate = odesolnext

  Uh = allocate_space(sol.trial)
  Uh = evaluate!(Uh, sol.trial, tF)
  uh = FEFunction(Uh, uF)

  state = (Uh, odesolstate)
  (uh, tF), state
end

function Base.iterate(sol::TransientFESolution, state)
  Uh, odesolstate = state

  odesolnext = iterate(sol.odesol, odesolstate)
  if isnothing(odesolnext)
    return nothing
  end

  (uF, tF), odesolstate = odesolnext

  Uh = evaluate!(Uh, sol.trial, tF)
  uh = FEFunction(Uh, uF)

  state = (Uh, odesolstate)
  (uh, tF), state
end

Base.IteratorSize(::Type{TransientFESolution}) = Base.SizeUnknown()

########
# Test #
########
"""
    test_transient_fe_solution(sol::TransientFESolution) -> Bool

Test the interface of `TransientFESolution` specializations.
"""
function test_transient_fe_solution(sol::TransientFESolution)
  for (uh_n, t_n) in sol
    @test t_n isa Real
    @test uh_n isa FEFunction
  end
  true
end

"""
    test_transient_fe_solver(
      solver::ODESolver, op::TransientFEOperator,
      u0, t0, tF
    ) -> Bool

Test the interface of `ODESolver` specializations on `TransientFEOperator`s.
"""
function test_transient_fe_solver(
  solver::ODESolver, op::TransientFEOperator,
  u0, t0, tF
)
  solution = solve(solver, op, u0, t0, tF)
  test_transient_fe_solution(solution)
end
