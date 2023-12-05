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
  odeop::ODEOperator{<:AbstractQuasilinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::ForwardEuler, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real
)
  ForwardEulerNonlinearOperator(
    odeop, odeopcache,
    t0, us0, dt
  )
end

function DiscreteODEOperator(
  odeslvr::ForwardEuler, odeop::ODEOperator{<:AbstractQuasilinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real
)
  J, r = disopcache
  ForwardEulerLinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    J, r
  )
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
    us = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us)
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

######################
# Nonlinear operator #
######################
"""
    struct ForwardEulerNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the forward Euler scheme:
```math
residual(t0, u0, v0) = 0,

t0 = t_n
u0 = u_n
v0 = x,

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

function Algebra.residual!(
  r::AbstractVector,
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t0, (u0,) = disop.t0, disop.us0
  us0 = (u0, x)
  residual!(r, disop.odeop, t0, us0, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t0, (u0,) = disop.t0, disop.us0
  us0 = (u0, x)
  allocate_jacobian(disop.odeop, t0, us0, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t0, (u0,) = disop.t0, disop.us0
  us0 = (u0, x)
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, disop.odeop, t0, us0, 1, 1, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ForwardEulerNonlinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, dt, us0 = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0)

  uF, = usF
  disslvrcache = solve!(uF, disslvr, disop, disslvrcache)

  # Express usF in terms of the solution of the discrete ODE operator
  usF = _finalize_euler!(usF, us0, uF, dt)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator

Linear discrete operator corresponding to the forward Euler scheme:
```math
residual(t0, u0, v0) = mass(t0, u0) v0 + res(t0, u0) = 0,

t0 = t_n
u0 = u_n
v0 = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::ForwardEulerLinearOperator) = disop.J

Algebra.get_vector(disop::ForwardEulerLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ForwardEulerLinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  t0, dt, us0 = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0)

  # Update jacobian and residual
  u0, = us0
  uF, = usF
  us = (u0, u0)
  w = 1
  filter = (true, true, false)

  fillstored!(J, zero(eltype(J)))
  jacobian!(J, odeop, t0, us, 1, w, odeopcache)
  residual!(r, odeop, t0, us, odeopcache; filter)
  rmul!(r, -1)

  # Solve the discrete ODE operator
  uF = usF[1]
  disslvrcache = solve!(uF, disslvr, disop, disslvrcache)

  # Express usF in terms of the solution of the discrete ODE operator
  usF = _finalize_euler!(usF, us0, uF, dt)

  (usF, disslvrcache)
end

#############
# Finalizer #
#############
function _finalize_euler!(
  usF::NTuple{1,AbstractVector}, us0::NTuple{1,AbstractVector},
  x::AbstractVector, dt::Real
)
  u0, = us0
  uF, = usF
  @. uF = u0 + dt * x
  (uF,)
end
