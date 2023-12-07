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
  t0::Real, us0::NTuple{2,AbstractVector}
)
  nothing
end

function allocate_disopcache(
  odeslvr::ForwardEuler,
  odeop::ODEOperator{<:AbstractQuasilinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::ForwardEuler, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
)
  ForwardEulerNonlinearOperator(
    odeop, odeopcache,
    t0, u0, dt
  )
end

function DiscreteODEOperator(
  odeslvr::ForwardEuler, odeop::ODEOperator{<:AbstractQuasilinearODE},
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
)
  J, r = disopcache
  ForwardEulerLinearOperator(
    odeop, odeopcache,
    t0, u0, dt,
    J, r
  )
end

function solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::ForwardEuler, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector},
  cache
)
  u0 = us0[1]
  dt = get_dt(odeslvr)
  tF = t0 + dt

  # Allocate or unpack cache
  if isnothing(cache)
    us0_full = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us0_full)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, us0_full)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create and solve discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop, odeopcache, disopcache,
    t0, u0, dt
  )

  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  cache = (odeopcache, disopcache, disslvrcache)

  (tF, usF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct ForwardEulerNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the forward Euler scheme:
```math
residual(tx, ux, vx) = 0,

tx = t_n
ux = u_n
vx = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ForwardEulerNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  u0::AbstractVector
  dt::Real
end

function Algebra.residual!(
  r::AbstractVector,
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t0, u0 = disop.t0, disop.u0
  # Residual: (u0, x)
  tx = t0
  usx = (u0, x)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t0, u0 = disop.t0, disop.u0
  tx = t0
  usx = (u0, x)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::ForwardEulerNonlinearOperator,
  x::AbstractVector
)
  t0, u0 = disop.t0, disop.u0
  # Jacobian: (u0, x)
  tx = t0
  usx = (u0, x)
  w = 1
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, disop.odeop, tx, usx, 1, w, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ForwardEulerNonlinearOperator,
  disslvrcache
)
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, t0)

  # Solve the discrete ODE operator
  x = uF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  uF = _finalize_euler!(uF, u0, dt)
  usF = (uF,)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator

Linear discrete operator corresponding to the forward Euler scheme:
```math
residual(tx, ux, vx) = mass(tx, ux) vx + res(tx, ux) = 0,

tx = t_n
ux = u_n
vx = x.
```
"""
struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  u0::AbstractVector
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
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt

  J, r = disop.J, disop.r

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, t0)

  # Jacobian: (u0, x)
  # Residual: (u0, x)
  # Trick: use uF to store 0
  x = uF
  fill!(x, zero(eltype(x)))

  tx = t0
  usx = (u0, x)
  w = 1
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, odeop, tx, usx, 1, w, odeopcache)
  residual!(r, odeop, tx, usx, odeopcache)
  rmul!(r, -1)

  # Solve the discrete ODE operator
  x = uF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  uF = _finalize_euler!(uF, u0, dt)
  usF = (uF,)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _finalize_euler!(uF::AbstractVector, u0::AbstractVector, dt::Real)
  # @. uF = u0 + dt * x
  axpby!(1, u0, dt, uF)
  uF
end
