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
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# Γ = BoundaryTriangulation(model, tags="boundary")
# dΓ = Measure(Γ, degree)
# nΓ = get_normal_vector(Γ)
# h = 1 / partition[1]
# β = 10

# ODE operator
m(t, u, v) = ∫(u ⋅ v) * dΩ
a(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ # - ∫(u ⋅ (nΓ ⋅ ∇(v)) + u ⋅ (nΓ ⋅ ∇(v)) - β / h * (u ⋅ v)) * dΓ
b(t, v) = ∫(f(t) ⋅ v) * dΩ # - ∫(u(t) ⋅ (nΓ ⋅ ∇(v)) - β / h * (u(t) ⋅ v)) * dΓ

mass(t, u, v) = m(t, ∂t(u), v)
stiffness(t, u, v) = a(t, u, v)
forcing(t, v) = (-1) * b(t, v)
res_quasilinear(t, u, v) = stiffness(t, u, v) + forcing(t, v)
res_nonlinear(t, u, v) = mass(t, u, v) + res_quasilinear(t, u, v)
jac(t, u, du, v) = a(t, du, v)
jac_t(t, u, dut, v) = m(t, dut, v)

feop_nonlinear = TransientFEOperator(res_nonlinear, jac, jac_t, U, V)
feop_quasilinear = TransientQuasilinearFEOperator(mass, res_quasilinear, jac, jac_t, U, V)
feop_semilinear = TransientSemilinearFEOperator(mass, res_quasilinear, jac, jac_t, U, V)
feop_linear = TransientLinearFEOperator(mass, stiffness, forcing, jac, jac_t, U, V)

feops = [
  feop_nonlinear,
  feop_quasilinear,
  feop_semilinear,
  feop_linear
]

# ODE solver
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)

disslvr = NLSolverMock()
odeslvr = ODESolverMock1(disslvr, dt)

for feop in feops
  fesltn = solve(odeslvr, feop, uh0, t0, tF)
  @test test_transient_fe_solution(fesltn)

  tol = 1.0e-6
  l = 0
  for (uh_n, t_n) in fesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(∫(eh_n * eh_n) * dΩ))
    @test e_n < tol

    l += 1
  end

  @test l == trunc(Int, ceil((tF - t0) / dt))
end

end # module TransientFESolutionsTests
