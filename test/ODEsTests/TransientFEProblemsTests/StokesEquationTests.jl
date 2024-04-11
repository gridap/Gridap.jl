module StokesEquationTests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
ut(t) = x -> VectorValue(x[1], x[2]) * t
u = TimeSpaceFunction(ut)

pt(t) = x -> (x[1] - x[2]) * t
p = TimeSpaceFunction(pt)

# Geometry
domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

# FE spaces
order = 2
reffe_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
V = FESpace(model, reffe_u, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)

reffe_p = ReferenceFE(lagrangian, Float64, order - 1)
Q = FESpace(model, reffe_p, conformity=:H1, constraint=:zeromean)
P = TrialFESpace(Q)

X = TransientMultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])

# Integration
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# FE operator
ft(t) = x -> ∂t(u)(t, x) - Δ(u)(t, x) + ∇(p)(t, x)
f = TimeSpaceFunction(ft)
g = ∇ ⋅ u
mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, (v, q)) = ∫(f(t) ⋅ v) * dΩ + ∫(g(t) * q) * dΩ

res(t, (u, p), (v, q)) = mass(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, (v, q)) - ∫(p * (∇ ⋅ v)) * dΩ + ∫((∇ ⋅ u) * q) * dΩ
jac(t, (u, p), (du, dp), (v, q)) = stiffness(t, du, v) - ∫(dp * (∇ ⋅ v)) * dΩ + ∫((∇ ⋅ du) * q) * dΩ
jac_t(t, (u, p), (dut, dpt), (v, q)) = mass(t, dut, v)

tfeop_nl_man = TransientFEOperator(res, (jac, jac_t), X, Y)
tfeop_nl_ad = TransientFEOperator(res, X, Y)
tfeops = (
  tfeop_nl_man,
  tfeop_nl_ad,
)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
P0 = P(t0)
ph0 = interpolate_everywhere(p(t0), P0)
X0 = X(t0)
xh0 = interpolate_everywhere([uh0, ph0], X0)
xhs0 = (xh0,)

# ODE Solver
tol = 1.0e-6
sysslvr_l = LUSolver()
sysslvr_nl = NLSolver(sysslvr_l, show_trace=false, method=:newton, iterations=10)
odeslvrs = (
  ThetaMethod(sysslvr_nl, dt, 0.5),
)

# Tests
for odeslvr in odeslvrs
  for tfeop in tfeops
    fesltn = solve(odeslvr, tfeop, t0, tF, xhs0)

    for (t_n, xhs_n) in fesltn
      eh_n = u(t_n) - xhs_n[1]
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol

      eh_n = p(t_n) - xhs_n[2]
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol
    end
  end
end

end # module StokesEquationTests
