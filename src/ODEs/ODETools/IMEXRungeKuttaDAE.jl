"""
Implicit-Explicit Runge-Kutta solver for Differential Algebraic Equations.

This struct defines a solver for the DAE system of the form

      M(u,t)du/dt = f(u,t) + g(u,u*,t),
      A(u,u*,t) = 0.

  where `f` is a nonlinear function of the differential variables `u` and `t`
  that will treated implicitly and `g` is a nonlinear function of `u`, the
  algebraic variables `u*` and `t` that will be treated explicitly. `A` is a
  nonlinear operator that defines the algebraic constraints.
  The DAE is solved using an implicit-explicit Runge-Kutta method.
"""
struct IMEXRungeKuttaDAE <: ODESolver
  nls_stage::NonlinearSolver
  nls_update::NonlinearSolver
  dae_solver::NonlinearSolver
  dt::Float64
  tableau::IMEXButcherTableau
  function IMEXRungeKuttaDAE(nls_stage::NonlinearSolver, nls_update::NonlinearSolver,
     dae_solver::NonlinearSolver, dt, type::Symbol)
    bt = IMEXButcherTableau(type)
    new(nls_stage, nls_update, dae_solver, dt, bt)
  end
end

"""
solve_step!(uf,odesol,op,dae_op,u0,t0,cache)

Solve one step of the DAE problem defined by `op` using the ODE solver `odesol`
  with initial solution `u0` at time `t0`. The solution is stored in `uf` and
  the final time in `tf`. The cache is used to store the solution of the
  nonlinear system of equations and auxiliar variables, including the algebraic
  variables.
"""
function solve_step!(uf::AbstractVector,
  solver::IMEXRungeKuttaDAE,
  op::DAEOperator,
  u0::AbstractVector,
  t0::Real,
  cache)

  # Unpack variables
  dt = solver.dt
  s = solver.tableau.s
  aᵢ = solver.tableau.aᵢ
  bᵢ = solver.tableau.bᵢ
  aₑ = solver.tableau.aₑ
  bₑ = solver.tableau.bₑ
  c = solver.tableau.c
  d = solver.tableau.d

  # Unpack operators
  ode_op = op.ode_op
  alg_op = op.alg_op

  # Create cache if not there
  if cache === nothing
    ode_cache = allocate_cache(ode_op)
    vi = similar(u0)
    fi = Vector{typeof(u0)}(undef,0)
    gi = Vector{typeof(u0)}(undef,0)
    u_alg = zero_initial_guess(alg_op)
    for i in 1:s
      push!(fi,similar(u0))
      push!(gi,similar(u0))
    end
    nls_stage_cache = nothing
    nls_update_cache = nothing
    dae_solver_cache = nothing
  else
    ode_cache, vi, fi, gi, u_alg, nls_stage_cache, nls_update_cache, dae_solver_cache = cache
  end

  # Create RKNL stage operator
  nlop_stage = IMEXRungeKuttaStageNonlinearOperator(ode_op,t0,dt,u0,ode_cache,vi,fi,gi,0,aᵢ,aₑ)

  # Compute intermediate stages
  for i in 1:s

    # Update time
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop_stage,ti,fi,gi,i)

    if(aᵢ[i,i]==0)
      # Skip stage solve if a_ii=0 => u_i=u_0, f_i = f_0, gi = g_0
      @. uf = u0
    else
      # solve at stage i
      nls_stage_cache = solve!(uf,solver.nls_stage,nlop_stage,nls_stage_cache)
    end

    # Solve algebraic constraints
    dae_solver_cache = solve!((u_alg,uf),solver.dae_solver,alg_op,dae_solver_cache)

    # Update RHS at stage i using solution at u_i
    rhs!(nlop_stage, uf)
    explicit_rhs!(nlop_stage, (u_alg,uf))

  end

  # Update final time
  tf = t0+dt

  # Skip final update if not necessary
  if !(c[s]==1.0 && aᵢ[s,:] == bᵢ && aₑ[s,:] == bₑ)

    # Create RKNL final update operator
    ode_cache = update_cache!(ode_cache,ode_op,tf)
    nlop_update = IMEXRungeKuttaUpdateNonlinearOperator(ode_op,tf,dt,u0,ode_cache,vi,fi,gi,s,bᵢ,bₑ)

    # solve at final update
    nls_update_cache = solve!(uf,solver.nls_update,nlop_update,nls_update_cache)

    # Solve algebraic constraints
    dae_solver_cache = solve!((u_alg,uf),solver.dae_solver,alg_op,dae_solver_cache)

  end

  # Update final cache
  cache = (ode_cache, vi, fi, gi, u_alg, nls_stage_cache, nls_update_cache, dae_solver_cache)

  return (uf, tf, cache)

end

function explicit_rhs!(op::RungeKuttaNonlinearOperator, x::Tuple{Vararg{AbstractVector}})
  u_alg, u = x
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  g = op.gi
  explicit_rhs!(g[op.i],op.odeop,op.ti,(u,v),u_alg,op.ode_cache)
end
