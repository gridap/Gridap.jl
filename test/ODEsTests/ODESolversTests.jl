module ODESolversTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

a, b, c = 1.0, 0.0, 1.0
ode_op = ODEOperatorMock{ConstantMassODE}(a, b, c, 1)

t0 = 0.0
tF = 1.0
dt = 0.1
u0 = 2 * ones(2)

α = 1 - dt * a
β = -dt * b
γ = 1 - dt * c
@assert !iszero(α)
@assert !iszero(γ)

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

D = α * γ
@test u[1] ≈ γ * u0[1] / D
@test u[2] ≈ (α * u0[2] - β * u0[1]) / D

# ODESolver tests
ode_solver = ODESolverMock(nls, dt)
uF = copy(u0)
fill!(uF, 1)

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

@test test_ode_solver(ode_solver, ode_op, u0, t0, tF)

end # module ODESolversTests
