using LinearAlgebra
using LinearAlgebra: fillstored!

using Gridap
using Gridap.Algebra
using Gridap.ODEs
using Gridap.Polynomials

#######################
# NonlinearSolverMock #
#######################
"""
    struct NonlinearSolverMock <: NonlinearSolver end

Mock `NonlinearSolver` for `NonlinearStageOperator` (simple Newton-Raphson with
Backslash) which defaults to Backslash for `linearStageOperator`.
"""
struct NonlinearSolverMock <: NonlinearSolver
  atol::Real
  rtol::Real
  maxiter::Integer
end

function Algebra.solve!(
  x::AbstractVector, nls::NonlinearSolverMock,
  nlop::NonlinearStageOperator, cache
)
  atol, rtol, maxiter = nls.atol, nls.rtol, nls.maxiter

  if isnothing(cache)
    J = allocate_jacobian(nlop, x)
    r = allocate_residual(nlop, x)
  else
    J, r = cache
  end

  # Update residual and check convergence
  residual!(r, nlop, x)
  nr_prev = norm(r)
  converged = (nr_prev < atol)

  iter = 0
  while !converged
    if iter > maxiter
      throw("NonlinearSolverMock did not converge.")
    end

    # Update x
    jacobian!(J, nlop, x)
    dx = J \ r
    axpy!(-1, dx, x)

    # Update residual and check convergence
    residual!(r, nlop, x)
    nr = norm(r)
    converged = (nr < atol) || (nr < rtol * nr_prev)
    nr_prev = nr

    iter += 1
  end

  (J, r)
end

function Algebra.solve!(
  x::AbstractVector, nls::NonlinearSolverMock,
  lop::LinearStageOperator, cache
)
  copy!(x, lop.J \ lop.r)
  rmul!(x, -1)
  cache
end

#################
# ODESolverMock #
#################
"""
    struct ODESolverMock <: ODESolver end

Mock `ODESolver` for ODEs of arbitrary order, using a backward Euler scheme.
```
res(tx, ux[0], ..., ux[N]) = 0,

tx = t_n + dt
ux[i] = ∂t^i[u](t_n) + ∑_{i + 1 ≤ j ≤ N} 1/(j - i)! dt^(j - i) ux[j].
ux[N] = x

∂t^i[u](t_(n+1)) = ux[i].
```
"""
struct ODESolverMock <: ODESolver
  sysslvr::NonlinearSolverMock
  dt::Real
end

##################
# Nonlinear case #
##################
function ODEs.allocate_odecache(
  odeslvr::ODESolverMock, odeop::ODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}
)
  u0 = us0[1]
  us0N = (us0..., u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  usx = copy.(us0)

  N = get_order(odeop)
  dt = odeslvr.dt
  ws = ntuple(i -> 0, N + 1)
  ws = Base.setindex(ws, 1, N + 1)
  for i in N:-1:1
    wi, coef = 0, 1
    for j in i+1:N+1
      coef = coef * dt / (j - i)
      wi += coef * ws[j]
    end
    ws = Base.setindex(ws, wi, i)
  end

  sysslvrcache = nothing
  odeslvrcache = (usx, ws, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ODEs.ode_march!(
  stateF::Tuple{Vararg{AbstractVector}},
  odeslvr::ODESolverMock, odeop::ODEOperator,
  t0::Real, state0::Tuple{Vararg{AbstractVector}},
  odecache
)
  # Unpack inputs
  us0, usF = state0, stateF
  odeslvrcache, odeopcache = odecache
  usx, ws, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt = odeslvr.dt

  # Define scheme
  tx = t0 + dt
  _usx(x) = _stage_mock!(usx, us0, dt, x)

  # Update ODE operator cache
  update_odeopcache!(odeopcache, odeop, tx)

  # Create and solve stage operator
  stageop = NonlinearStageOperator(odeop, odeopcache, tx, _usx, ws)

  x = usF[1]
  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  tF = t0 + dt
  copy!(usx[1], x)
  x = usx[1]
  usF = _convert_mock!(usF, us0, dt, x)
  stateF = usF

  # Pack outputs
  odeslvrcache = (usx, ws, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _stage_mock!(
  usx::Tuple{Vararg{AbstractVector}}, us0::Tuple{Vararg{AbstractVector}},
  dt::Real, x::AbstractVector
)
  _convert_mock!(usx, us0, dt, x)
  (usx..., x)
end

function _convert_mock!(
  usF::Tuple{Vararg{AbstractVector}}, us0::Tuple{Vararg{AbstractVector}},
  dt::Real, x::AbstractVector
)
  # usF[i] = us0[i] + ∑_{i < j ≤ N+1} 1/(j - i)! dt^(j - i) usF[j]
  N = length(us0)
  for i in N:-1:1
    ui0, uiF = us0[i], usF[i]
    copy!(uiF, ui0)
    coef = 1
    for j in i+1:N
      coef = coef * dt / (j - i)
      axpy!(coef, usF[j], uiF)
    end
    j = N + 1
    coef = coef * dt / (j - i)
    axpy!(coef, x, uiF)
  end
  usF
end
