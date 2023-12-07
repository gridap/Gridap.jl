using LinearAlgebra
using LinearAlgebra: fillstored!

using Gridap
using Gridap.Algebra
using Gridap.ODEs
using Gridap.Polynomials

#########################
# DiscreteODESolverMock #
#########################
"""
    struct DiscreteODESolverMock <: NonlinearSolver end

Mock `DiscreteODESolver` (simple Newton-Raphson with Backslash).
"""
struct DiscreteODESolverMock <: NonlinearSolver
  atol::Real
  rtol::Real
  maxiter::Integer
end

function Algebra.solve!(
  x::AbstractVector, disslvr::DiscreteODESolverMock,
  op::NonlinearOperator, cache
)
  atol, rtol, maxiter = disslvr.atol, disslvr.rtol, disslvr.maxiter

  if isnothing(cache)
    r = allocate_residual(op, x)
    J = allocate_jacobian(op, x)
  else
    r, J = cache
  end

  # Update residual and check convergence
  residual!(r, op, x)
  nr_prev = norm(r)
  converged = (nr_prev < atol)

  iter = 0
  while !converged
    if iter > maxiter
      throw("DiscreteODESolverMock did not converge.")
    end

    # Update jacobian
    jacobian!(J, op, x)

    # Update x
    dx = J \ r
    axpy!(-1, dx, x)

    # Update residual and check convergence
    residual!(r, op, x)
    nr = norm(r)
    converged = (nr < atol) || (nr < rtol * nr_prev)
    nr_prev = nr

    iter += 1
  end

  (r, J)
end

#################
# ODESolverMock #
#################
"""
    struct ODESolverMock <: ODESolver end

Mock `ODESolver` for ODEs of arbitrary order, using a backward Euler scheme.
"""
struct ODESolverMock <: ODESolver
  disslvr::DiscreteODESolverMock
  dt::Real
end

function ODEs.get_dt(odeslvr::ODESolverMock)
  odeslvr.dt
end

function ODEs.allocate_disopcache(
  odeslvr::ODESolverMock,
  odeop::ODEOperator, odeopcache,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}
)
  N = get_order(odeop)
  ntuple(k -> zero(us0[k]), N)
end

function ODEs.DiscreteODEOperator(
  odeslvr::ODESolverMock, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}, dt::Real, tx::Real
)
  usx = disopcache

  N = get_order(odeop)
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

  DiscreteODEOperatorMock(
    odeop, odeopcache,
    us0, dt, tx, usx, ws
  )
end

function ODEs.solve_step!(
  usF::Tuple{Vararg{AbstractVector}},
  odeslvr::ODESolverMock, odeop::ODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}},
  cache
)
  dt = get_dt(odeslvr)
  tF = t0 + dt

  # Allocate or unpack cache
  if isnothing(cache)
    us0_full = (us0..., us0[1])
    odeopcache = allocate_odeopcache(odeop, t0, us0_full)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, us0_full)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create and solve discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop, odeopcache, disopcache,
    t0, us0, dt, tF
  )

  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  cache = (odeopcache, disopcache, disslvrcache)

  (tF, usF, cache)
end

###########################
# DiscreteODEOperatorMock #
###########################
"""
    struct DiscreteODEOperatorMock <: DiscreteODEOperator end

Mock backward Euler operator for arbitrary-order ODEs:
```math
res(tx, ux[0], ..., ux[N]) = 0,

tx = t_n + dt
ux[i] = ∂t^i[u](t_n) + ∑_{i + 1 ≤ j ≤ N} 1/(j - i)! dt^(j - i) ux[j].
ux[N] = x

∂t^i[u](t_(n+1)) = ux[i].
```
"""
struct DiscreteODEOperatorMock <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  us0::Tuple{Vararg{AbstractVector}}
  dt::Real
  tx::Real
  usx::Tuple{Vararg{AbstractVector}}
  ws::Tuple{Vararg{Real}}
end

function Algebra.allocate_residual(
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  dt, tx = disop.dt, disop.tx
  us0, usx = disop.us0, disop.usx
  usx = _stage_mock!(usx, dt, us0, x)
  allocate_residual(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, dt = disop.us0, disop.dt
  tx, usx = disop.tx, disop.usx
  usx = _stage_mock!(usx, dt, us0, x)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
  r
end

function Algebra.allocate_jacobian(
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, dt = disop.us0, disop.dt
  tx, usx = disop.tx, disop.usx
  usx = _stage_mock!(usx, dt, us0, x)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, dt = disop.us0, disop.dt
  tx, usx = disop.tx, disop.usx
  usx = _stage_mock!(usx, dt, us0, x)
  ws = disop.ws
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
end

function Algebra.solve!(
  usF::Tuple{Vararg{AbstractVector}},
  disslvr::NonlinearSolver, disop::DiscreteODEOperatorMock,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  us0, dt = disop.us0, disop.dt
  tx = disop.tx
  x = usF[1]

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tx)

  # Solve the discrete ODE operator
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  copy!(disop.usx[1], x)
  x = disop.usx[1]
  usF = _finalize_mock!(usF, dt, us0, x)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _stage_mock!(
  usF::Tuple{Vararg{AbstractVector}}, dt::Real,
  us0::Tuple{Vararg{AbstractVector}}, x::AbstractVector
)
  usF = _finalize_mock!(usF, dt, us0, x)
  usF = (usF..., x)
  usF
end

function _finalize_mock!(
  usF::Tuple{Vararg{AbstractVector}}, dt::Real,
  us0::Tuple{Vararg{AbstractVector}}, x::AbstractVector
)
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
