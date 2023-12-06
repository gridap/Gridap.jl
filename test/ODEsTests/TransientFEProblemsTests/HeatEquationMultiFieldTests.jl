module HeatEquationMultifieldTests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * (1 + t)
u(t::Real) = x -> u(x, t)
u(x) = t -> u(x, t)

∂tu(x, t) = ∂t(u)(x, t)
∂tu(t::Real) = x -> ∂tu(x, t)

# Geometry
domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

# FE spaces
order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V = FESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)

Y = MultiFieldFESpace([V, V])
X = TransientMultiFieldFESpace([U, U])

# Integration
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# FE operator
f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)
mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

_res(t, u, v) = mass(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
res(t, (u1, u2), (v1, v2)) = _res(t, u1, v1) + _res(t, u2, v2)
jac(t, us, (du1, du2), (v1, v2)) = stiffness(t, du1, v1) + stiffness(t, du2, v2)
jac_t(t, us, (dut1, dut2), (v1, v2)) = mass(t, dut1, v1) + mass(t, dut2, v2)

feop = TransientFEOperator(res, jac, jac_t, X, Y)
feop_AD = TransientFEOperator(res, X, Y)
feops = (feop, feop_AD,)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
X0 = X(t0)
xh0 = interpolate_everywhere([uh0, uh0], X0)
xhs0 = (xh0,)

# ODE Solver
tol = 1.0e-6
disslvr = LUSolver()
odeslvrs = (
  MidPoint(disslvr, dt),
)

# Tests
for odeslvr in odeslvrs
  for feop in feops
    fesltn = solve(odeslvr, feop, xhs0, t0, tF)

    for (xhs_n, t_n) in fesltn
      eh_n = u(t_n) - xhs_n[1]
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol

      eh_n = u(t_n) - xhs_n[2]
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol
    end
  end
end

end # module HeatEquationMultifieldTests
