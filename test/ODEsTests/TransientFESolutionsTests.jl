module TransientFESolutionsTests

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
partition = (2, 2)
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

# Γ = BoundaryTriangulation(model, tags="boundary")
# dΓ = Measure(Γ, degree)
# nΓ = get_normal_vector(Γ)
# h = 1 / partition[1]
# β = 10

# ODE operator
m(t, u, v) = ∫(u * v) * dΩ
a(t, u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ # - ∫(v ⋅ (nΓ ⋅ ∇(u)) + u ⋅ (nΓ ⋅ ∇(v)) - β / h * (v ⋅ u)) * dΓ
b(t, v) = ∫(v * f(t)) * dΩ # - ∫(u(t) ⋅ (nΓ ⋅ ∇(v)) - β / h * (v ⋅ u(t))) * dΓ

res(t, u, v) = m(t, ∂t(u), v) + a(t, u, v) - b(t, v)
jac(t, u, du, v) = a(t, du, v)
jac_t(t, u, dut, v) = m(t, dut, v)
op = TransientFEOperator(res, jac, jac_t, U, V)

# ODE solver
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)

nls = NLSolverMock()
ode_solver = ODESolverMock(nls, dt)

sol = solve(ode_solver, op, uh0, t0, tF)
@test test_transient_fe_solution(sol)

tol = 1.0e-6
l = 0
for (uh_n, t_n) in sol
  eh_n = u(t_n) - uh_n
  e_n = sqrt(sum(∫(eh_n * eh_n) * dΩ))
  @test e_n < tol

  global l
  l += 1
end

@test l == trunc(Int, ceil((tF - t0) / dt))

end # module TransientFESolutionsTests
