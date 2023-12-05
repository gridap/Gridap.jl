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
  t::Real, x::AbstractVector
)
  N = get_order(odeop)
  ntuple(i -> zero(x), N + 1)
end

function ODEs.DiscreteODEOperator(
  odeslvr::ODESolverMock, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}, dt::Real
)
  tx = t0 + dt
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
    us0, usx, tx, dt, ws
  )
end

function ODEs.solve_step!(
  usF::Tuple{Vararg{AbstractVector}},
  odeslvr::ODESolverMock, odeop::ODEOperator,
  us0::Tuple{Vararg{AbstractVector}}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0 = first(us0)
  dt = get_dt(odeslvr)

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, (us0..., u0))
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    t0, us0, dt
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
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
ux[N] = x
ux[i] = ∂t^i[u](t_n) + ∑_{i + 1 ≤ j ≤ N} 1/(j - i)! dt^(j - i) ux[j],

∂t^i[u](t_(n+1)) = ux[i].
```
"""
struct DiscreteODEOperatorMock <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  us0::Tuple{Vararg{AbstractVector}}
  usx::Tuple{Vararg{AbstractVector}}
  tx::Real
  dt::Real
  ws::Tuple{Vararg{Real}}
end

function Algebra.allocate_residual(
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, usx = disop.us0, disop.usx
  tx, dt = disop.tx, disop.dt
  usx = _set_usx!(usx, us0, x, dt)
  allocate_residual(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, usx = disop.us0, disop.usx
  tx, dt = disop.tx, disop.dt
  usx = _set_usx!(usx, us0, x, dt)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, usx = disop.us0, disop.usx
  tx, dt = disop.tx, disop.dt
  usx = _set_usx!(usx, us0, x, dt)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::DiscreteODEOperatorMock,
  x::AbstractVector
)
  us0, usx = disop.us0, disop.usx
  tx, dt = disop.tx, disop.dt
  usx = _set_usx!(usx, us0, x, dt)
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
  us0, usx = disop.us0, disop.usx
  tx, dt = disop.tx, disop.dt

  # Update dirichlet boundary conditions
  update_odeopcache!(odeopcache, odeop, tx)

  # Solve discrete ODE operator
  x = last(usF)
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalizer
  _set_usx!(usx, us0, x, dt)
  for i in eachindex(usF)
    uiF, uix = usF[i], usx[i]
    copy!(uiF, uix)
    usF = Base.setindex(usF, uiF, i)
  end

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _set_usx!(usx, us0, x, dt)
  N = length(us0)

  uNx = usx[N+1]
  copy!(uNx, x)
  usx = Base.setindex(usx, uNx, N + 1)

  for i in N:-1:1
    ui0, uix = us0[i], usx[i]
    copy!(uix, ui0)
    coef = 1
    for j in i+1:N+1
      coef = coef * dt / (j - i)
      axpy!(coef, usx[j], uix)
    end
    usx = Base.setindex(usx, uix, i)
  end

  usx
end
