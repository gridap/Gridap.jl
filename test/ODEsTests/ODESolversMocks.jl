using LinearAlgebra
using LinearAlgebra: fillstored!

using Gridap
using Gridap.Algebra
using Gridap.ODEs
using Gridap.ODEs: _u_from_v!

################
# NLSolverMock #
################
"""
    struct NLSolverMock <: NonlinearSolver end

Mock linear system solver.
"""
struct NLSolverMock <: NonlinearSolver
end

function Algebra.solve!(
  x::AbstractVector, nls::NLSolverMock,
  op::NonlinearOperator, cache::Nothing
)
  r = residual(op, x)
  J = jacobian(op, x)

  dx = -J \ r
  axpy!(1, dx, x)

  cache = (r, J, dx)
  cache
end

function Algebra.solve!(
  x::AbstractVector, nls::NLSolverMock,
  op::NonlinearOperator, cache
)
  r, J, dx = cache

  residual!(r, op, x)
  jacobian!(J, op, x)

  dx = -J \ r
  axpy!(1, dx, x)

  cache
end

##################
# ODESolverMock1 #
##################
"""
    struct ODESolverMock1 <: ODESolver end

Mock first-order ODE solver (backward Euler).
"""
struct ODESolverMock1 <: ODESolver
  disslvr::NLSolverMock
  dt::Float64
end

function ODEs.get_dt(odeslvr::ODESolverMock1)
  odeslvr.dt
end

function ODEs.allocate_disopcache(
  odeslvr::ODESolverMock1,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  nothing
end

function ODEs.DiscreteODEOperator(
  odeslvr::ODESolverMock1, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real
)
  DiscreteODEOperatorMock1(odeop, odeopcache, t0, us0, dt)
end

function ODEs.solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::ODESolverMock1, odeop::ODEOperator,
  us0::NTuple{1,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, = us0
  dt = get_dt(odeslvr)

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, (u0, u0))
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

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
end

############################
# DiscreteODEOperatorMock1 #
############################
"""
    struct DiscreteODEOperatorMock1 <: DiscreteODEOperator end

Mock backward Euler operator for first-order ODEs:
```math
res(t_n + dt, u_n + dt * x, x) = 0.
u_n+1 = u_n + dt * x
```
"""
struct DiscreteODEOperatorMock1 <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
end

function Algebra.allocate_residual(
  disop::DiscreteODEOperatorMock1,
  x::AbstractVector
)
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* x
  allocate_residual(disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::DiscreteODEOperatorMock1,
  x::AbstractVector
)
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* x
  residual!(r, disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::DiscreteODEOperatorMock1,
  x::AbstractVector
)
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* x
  allocate_jacobian(disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::DiscreteODEOperatorMock1,
  x::AbstractVector
)
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* x
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, t, (u, x), (disop.dt, 1), disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::DiscreteODEOperatorMock1,
  disslvrcache
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0 + dt)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  _u_from_v!(uF, u0, dt, vF)

  usF = (uF,)
  (usF, disslvrcache)
end

##################
# ODESolverMock2 #
##################
"""
    struct ODESolverMock2 <: ODESolver end

Mock second-order ODE solver (backward Euler).
"""
struct ODESolverMock2 <: ODESolver
  disslvr::NLSolverMock
  dt::Float64
end

function ODEs.get_dt(odeslvr::ODESolverMock2)
  odeslvr.dt
end

function ODEs.allocate_disopcache(
  odeslvr::ODESolverMock2,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  nothing
end

function ODEs.DiscreteODEOperator(
  odeslvr::ODESolverMock2, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real
)
  DiscreteODEOperatorMock2(odeop, odeopcache, t0, us0, dt)
end

function ODEs.solve_step!(
  usF::NTuple{2,AbstractVector},
  odeslvr::ODESolverMock2, odeop::ODEOperator,
  us0::NTuple{2,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, _ = us0
  dt = get_dt(odeslvr)

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, (u0, u0, u0))
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

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
end

############################
# DiscreteODEOperatorMock2 #
############################
"""
    struct DiscreteODEOperatorMock2 <: DiscreteODEOperator end

Mock backward Euler operator for second-order ODEs:
```math
res(t_n + dt, u_n + dt * v_n + 3 * dt^2 / 2 * x, v_n + dt * x, x) = 0,
u_n+1 = u_n + dt * v_n + 3 * dt^2 / 2 * x
v_n+1 = v_n + dt * x
```
"""
struct DiscreteODEOperatorMock2 <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{2,AbstractVector}
  dt::Real
end

function Algebra.allocate_residual(
  disop::DiscreteODEOperatorMock2,
  x::AbstractVector
)
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* v0 .+ 3 * dt^2 / 2 .* x
  v = v0 .+ dt .* x
  allocate_residual(disop.odeop, t, (u, v, x), disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::DiscreteODEOperatorMock2,
  x::AbstractVector
)
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* v0 .+ 3 * dt^2 / 2 .* x
  v = v0 .+ dt .* x
  residual!(r, disop.odeop, t, (u, v, x), disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::DiscreteODEOperatorMock2,
  x::AbstractVector
)
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* v0 .+ 3 * dt^2 / 2 .* x
  v = v0 .+ dt .* x
  allocate_jacobian(disop.odeop, t, (u, v, x), disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::DiscreteODEOperatorMock2,
  x::AbstractVector
)
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0
  t = t0 + dt
  u = u0 .+ dt .* v0 .+ 3 * dt^2 / 2 .* x
  v = v0 .+ dt .* x
  fillstored!(J, zero(eltype(J)))
  ws = (3 * dt^2 / 2, disop.dt, 1)
  jacobians!(J, disop.odeop, t, (u, v, x), ws, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::DiscreteODEOperatorMock2,
  disslvrcache
)
  uF, vF = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0 + dt)

  aF = vF
  fill!(aF, zero(eltype(aF)))
  disslvrcache = solve!(aF, disslvr, disop, disslvrcache)

  @. uF = u0 + dt * v0 + 3 * dt^2 / 2 * aF
  @. vF = v0 + dt * aF

  usF = (uF, vF)
  (usF, disslvrcache)
end
