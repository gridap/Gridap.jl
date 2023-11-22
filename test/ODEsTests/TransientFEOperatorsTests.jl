module TransientFEOperatorsTests

using Test

using LinearAlgebra
using ForwardDiff

using Gridap
using Gridap.FESpaces
using Gridap.ODEs

include("ODESolversMocks.jl")

# Analytical functions
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * (t + 3.0)
u(t::Real) = x -> u(x, t)
u(x) = t -> u(x, t)

∂tu(x, t) = ∂t(u)(x, t)
∂tu(t::Real) = x -> ∂tu(x, t)

f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)

# Geometry
domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

# FE spaces
order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V = FESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)

# Integration
degree = 2 * order

Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# ODE operator
m(t, u, v) = ∫(u * v) * dΩ
a(t, u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ
b(t, v) = ∫(v * f(t)) * dΩ

res(t, u, v) = m(t, ∂t(u), v) + a(t, u, v) - b(t, v)
jac(t, u, du, v) = a(t, du, v)
jac_t(t, u, dut, v) = m(t, dut, v)
op = TransientFEOperator(res, jac, jac_t, U, V)

t0 = 0.0
t1 = 1.0
dt = 0.1

U0 = U(t0)
@test test_transient_fe_operator(op, zero(U0))

uh = interpolate_everywhere(t0, U0)

# FE operator at t=t0
dt⁻¹ = inv(dt)
_res(u, v) = dt⁻¹ * m(t0, u, v) + a(t0, u, v) - b(t0, v)
_jac(u, du, v) = dt⁻¹ * m(t0, du, v) + a(t0, du, v)
_op = FEOperator(_res, _jac, U0, V)

# Evaluate the residual and jacobian with TransientFEOperator
ode_op = get_algebraic_operator(op)
ode_cache = allocate_cache(ode_op)
uₜ = TransientCellField(uh, (uh,))

r = allocate_residual(op, t0, uh, ode_cache)
J = allocate_jacobian(op, t0, uh, ode_cache)
residual!(r, op, t0, uₜ, ode_cache)
jacobian!(J, op, t1, uₜ, 0, 1, ode_cache)
jacobian!(J, op, t1, uₜ, 1, dt⁻¹, ode_cache)

# Evaluate the residual and jacobian with FEOperator
_r = residual(_op, uh)
_J = jacobian(_op, uh)

@test all(r .≈ _r)
@test all(J .≈ _J)

end # module TransientFEOperatorsTests
