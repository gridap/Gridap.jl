"""
    struct ForwardEuler <: ODESolver end

Forward Euler ODE solver.
"""
struct ForwardEuler <: ODESolver
  disslvr::NonlinearSolver
  dt::Real
end

# ODESolver interface
function get_dt(odeslvr::ForwardEuler)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::ForwardEuler,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  nothing
end

function allocate_disopcache(
  odeslvr::ForwardEuler,
  odeop::ODEOperator{<:AbstractMassLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  J = allocate_jacobian(odeop, t, (x, x), odeopcache)
  r = allocate_residual(odeop, t, (x, x), odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::ForwardEuler, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real
)
  ForwardEulerNonlinearOperator(odeop, odeopcache, t0, us0, dt)
end

function DiscreteODEOperator(
  odeslvr::ForwardEuler, odeop::ODEOperator{<:AbstractMassLinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real
)
  J, r = disopcache
  ForwardEulerLinearOperator(odeop, odeopcache, J, r, t0, us0, dt)
end

function solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::ForwardEuler, odeop::ODEOperator,
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

######################
# Nonlinear operator #
######################
"""
    struct ForwardEulerNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the forward Euler scheme:
```math
residual(t_n, u_n, v_n) = 0,

v_n = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ForwardEulerNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
end

function Algebra.allocate_residual(
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t, (u,) = disop.t0, disop.us0
  allocate_residual(disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t, (u,) = disop.t0, disop.us0
  residual!(r, disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t, (u,) = disop.t0, disop.us0
  allocate_jacobian(disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t, (u,) = disop.t0, disop.us0
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, disop.odeop, t, (u, x), 1, 1, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ForwardEulerNonlinearOperator,
  disslvrcache
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  _u_from_v!(uF, u0, dt, vF)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator

Linear discrete operator corresponding to the forward Euler scheme:
```math
residual(t_n, u_n, v_n) = mass(t_n, u_n) v_n + res(t_n, u_n) = 0,

v_n = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
end

Algebra.get_matrix(disop::ForwardEulerLinearOperator) = disop.J

Algebra.get_vector(disop::ForwardEulerLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ForwardEulerLinearOperator,
  disslvrcache
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0)

  u, x = u0, u0
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, odeop, t0, (u, x), 1, 1, odeopcache)
  residual!(r, odeop, t0, (u, x), odeopcache, include_mass=false)
  rmul!(r, -1)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  _u_from_v!(uF, u0, dt, vF)

  (usF, disslvrcache)
end
