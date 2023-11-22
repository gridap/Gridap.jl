module ForwardEulerTests

using Test

using LinearAlgebra
using ForwardDiff

using Gridap
using Gridap.FESpaces
using Gridap.FESpaces: get_algebraic_operator
using Gridap.ODEs

import Gridap.ODEs: ∂t

# Analytical functions
# u(x, t) = (x[1] + x[2]) * t
# u(x, t) = (2 * x[1] + x[2]) * t
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * t
u(t::Real) = x -> u(x, t)
v(x) = t -> u(x, t)

∂tu(t) = x -> ForwardDiff.derivative(v(x), t)
∂tu(x, t) = ∂tu(t)(x)
∂t(::typeof(u)) = ∂tu

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

Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# ODE operator
res(t, u, v) = ∫(∂t(u) ⋅ v + ∇(u) ⊙ ∇(v) - f(t) ⋅ v) * dΩ
jac(t, u, du, v) = ∫(∇(du) ⊙ ∇(v)) * dΩ
jac_t(t, u, dut, v) = ∫(dut ⋅ v) * dΩ
op = TransientFEOperator(res, jac, jac_t, U, V)

# ODE solver
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)

function test_solver(ode_solver, tol)
  sol = solve(ode_solver, op, uh0, t0, tF)

  for (uh_n, t_n) in sol
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(∫(eh_n * eh_n) * dΩ))
    @test e_n < tol
  end
end

tol = 1.0e-4
ls = LUSolver()

ode_solvers = [
  ForwardEuler(ls, dt),
  ThetaMethod(ls, dt, 0.2),
  ThetaMethod(ls, dt, 0.5),
  ThetaMethod(ls, dt, 0.8),
  BackwardEuler(ls, dt)
]
for ode_solver in ode_solvers
  test_solver(ode_solver, tol)
end

end # module ForwardEulerTests
