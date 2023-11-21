module ODESolversTests

using Test

using Gridap
using Gridap.ODEs
using Gridap.ODEs.ODETools

include("ODEOperatorMocks.jl")
include("ODESolverMocks.jl")

a, b, c = 1.0, 0.0, 1.0
ode_op = ODEOperatorMock{ConstantMassODE}(a, b, c, 1)

t0 = 0.0
tF = 1.0
dt = 0.1
u0 = 2 * ones(2)

# OperatorMock tests
nl_op = OperatorMock(ode_op, u0, dt, tF, nothing)

u = zero_initial_guess(nl_op)
fill!(u, 1)
u̇ = (u - u0) ./ dt
dt⁻¹ = inv(dt)

r = allocate_residual(nl_op, u)
J = allocate_jacobian(nl_op, u)
residual!(r, nl_op, u)
jacobian!(J, nl_op, u)

@test r[1] ≈ u̇[1] - a
@test r[2] ≈ u̇[2] - b - c
@test J[1, 1] ≈ dt⁻¹ - a
@test J[1, 2] ≈ 0
@test J[2, 1] ≈ -b
@test J[2, 2] ≈ dt⁻¹ - c

# NLSolverMock tests
nls = NLSolverMock()

nls_cache = solve!(u, nls, nl_op)

r, J, du = nls_cache
@test r[1] ≈ u̇[1] - a
@test r[2] ≈ u̇[2] - b - c
@test J[1, 1] ≈ dt⁻¹ - a
@test J[1, 2] ≈ 0
@test J[2, 1] ≈ -b
@test J[2, 2] ≈ dt⁻¹ - c

D = (1 - dt * a) * (1 - dt * c)
@test u[1] ≈ (1 - dt * c) * u0[1] / D
@test u[2] ≈ ((dt * b) * u0[1] + (1 - dt * a) * u0[2]) / D

# ODESolver tests
ode_solver = ODESolverMock(nls, dt)
uF = copy(u0)
@. uF = 1

uF, tF, cache = solve_step!(uF, ode_solver, ode_op, u0, t0, nothing)

@test tF ≈ t0 + dt
@test all(uF .≈ u)

# ODESolutions tests
tF = 10.0
sol = GenericODESolution(ode_solver, ode_op, u0, t0, tF)
utF, state = Base.iterate(sol)
uF, tF = utF

@test tF ≈ t0 + dt
@test all(uF .≈ u)

end # module module ODESolversTests
