module ODESolversTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

t0 = randn()
tF = t0 + rand()
dt = (tF - t0) / 10
u0 = randn(2)

M = randn(2, 2)
K = randn(2, 2)
while iszero(det(M + dt * K))
  M = randn(2, 2)
  K = randn(2, 2)
end
f(t) = [cospi(t), sinpi(t)]
ode_op = ODEOperatorMock{MassLinearODE}(M, K, f)

# MockDiscreteODEOperator tests
dop = MockDiscreteODEOperator(ode_op, nothing, t0, u0, dt)

v = randn(2)

r = allocate_residual(dop, v)
J = allocate_jacobian(dop, v)
residual!(r, dop, v)
jacobian!(J, dop, v)

_r = M * v + K * (u0 + dt * v) + f(t0 + dt)
_J = M + dt * K
@test r ≈ _r
@test J ≈ _J

# NLSolverMock tests
sol = NLSolverMock()
sol_cache = solve!(v, sol, dop)

r, J, du = sol_cache
@test r ≈ _r
@test J ≈ _J

_v = -_J \ (K * u0 + f(t0 + dt))
@test v ≈ _v

# ODESolver tests
ode_solver = ODESolverMock(sol, dt)
uF = copy(u0)
uF, tF, cache = solve_step!(uF, ode_solver, ode_op, u0, t0, nothing)

@test tF ≈ t0 + dt
@test uF ≈ u0 + dt * _v

@test test_ode_solver(ode_solver, ode_op, t0, u0, t0, u0, dt)

end # module ODESolversTests
