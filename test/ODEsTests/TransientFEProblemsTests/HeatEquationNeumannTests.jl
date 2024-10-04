module HeatEquationNeumannTests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
ut(t) = x -> x[1] * (1 - x[2]) * (1 + t)
u = TimeSpaceFunction(ut)

# Geometry
domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)
dirichlet_tags = [1, 2, 3, 4, 5, 6]
neumanntags = [7, 8]

# FE spaces
order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V = FESpace(model, reffe, conformity=:H1, dirichlet_tags=dirichlet_tags)
U = TransientTrialFESpace(V, u)

# Integration
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

Γ = BoundaryTriangulation(model, tags=neumanntags)
dΓ = Measure(Γ, degree)
nΓ = get_normal_vector(Γ)

# FE operator
ft(t) = x -> ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)
mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
mass(t, u, ∂ₜu, v) = mass(t, ∂ₜu, v)
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing_Ω(t, v) = ∫(f(t) ⋅ v) * dΩ
forcing_Γ(t, v) = ∫((∇(u(t)) ⋅ nΓ) ⋅ v) * dΓ
forcing(t, v) = forcing_Ω(t, v) + forcing_Γ(t, v)

res(t, u, v) = mass(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = mass(t, dut, v)

res_ql(t, u, v) = stiffness(t, u, v) - forcing(t, v)
res_l(t, v) = forcing(t, v)

args_man = ((jac, jac_t), U, V)
tfeop_nl_man = TransientFEOperator(res, args_man...)
tfeop_ql_man = TransientQuasilinearFEOperator(mass, res_ql, args_man...)
tfeop_sl_man = TransientSemilinearFEOperator(mass, res_ql, args_man...)
tfeop_l_man = TransientLinearFEOperator((stiffness, mass), res_l, args_man...)

args_ad = (U, V)
tfeop_nl_ad = TransientFEOperator(res, args_ad...)
tfeop_ql_ad = TransientQuasilinearFEOperator(mass, res_ql, args_ad...)
tfeop_sl_ad = TransientSemilinearFEOperator(mass, res_ql, args_ad...)
tfeop_l_ad = TransientLinearFEOperator((stiffness, mass), res_l, args_ad...)

tfeops = (
  tfeop_nl_man,
  tfeop_ql_man,
  tfeop_sl_man,
  tfeop_l_man,
  tfeop_nl_ad,
  tfeop_ql_ad,
  tfeop_sl_ad,
  tfeop_l_ad
)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
uhs0 = (uh0,)

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
    fesltn = solve(odeslvr, tfeop, t0, tF, uhs0)

    for (t_n, uh_n) in fesltn
      eh_n = u(t_n) - uh_n
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol
    end
  end
end

end # module HeatEquationNeumannTests
