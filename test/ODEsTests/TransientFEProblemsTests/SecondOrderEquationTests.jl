module Order2FETests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs
# Analytical functions
ut(t) = x -> (1 - x[1]) * x[2] * (t^2 + 3.0)
u = TimeSpaceFunction(ut)

# Geometry
domain = (0, 1, 0, 1)
partition = (5, 5,)
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

# FE operator
order = 2
ft(t) = x -> ∂tt(u)(t, x) + ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)
mass(t, ∂ₜₜu, v) = ∫(∂ₜₜu ⋅ v) * dΩ
mass(t, u, ∂ₜₜu, v) = mass(t, ∂ₜₜu, v)
damping(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

res(t, u, v) = mass(t, u, ∂tt(u), v) + damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = damping(t, dut, v)
jac_tt(t, u, dutt, v) = mass(t, dutt, v)

res_ql(t, u, v) = damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
res_l(t, v) = forcing(t, v)

res0(t, u, v) = ∫(0 * u * v) * dΩ
jac0(t, u, du, v) = ∫(0 * du * v) * dΩ

args_man = ((jac, jac_t, jac_tt), U, V)
args0_man = ((jac0, jac0), U, V)
tfeop_nl_man = TransientFEOperator(res, args_man...)
tfeop_ql_man = TransientQuasilinearFEOperator(mass, res_ql, args_man...)
tfeop_sl_man = TransientSemilinearFEOperator(mass, res_ql, args_man...)
tfeop_l_man = TransientLinearFEOperator((stiffness, damping, mass), res_l, args_man...)

tfeop_im_man = TransientSemilinearFEOperator(mass, res_ql, args_man...)
tfeop_ex_man = TransientFEOperator(res0, args0_man...)
tfeop_imex_man = TransientIMEXFEOperator(tfeop_im_man, tfeop_ex_man)

args_ad = (U, V)
tfeop_nl_ad = TransientFEOperator(res, args_ad..., order=2)
tfeop_ql_ad = TransientQuasilinearFEOperator(mass, res_ql, args_ad..., order=2)
tfeop_sl_ad = TransientSemilinearFEOperator(mass, res_ql, args_ad..., order=2)
tfeop_l_ad = TransientLinearFEOperator((stiffness, damping, mass), res_l, args_ad...)

tfeop_im_ad = TransientSemilinearFEOperator(mass, res_ql, args_ad..., order=2)
tfeop_ex_ad = TransientFEOperator(res0, args_ad..., order=1)
tfeop_imex_ad = TransientIMEXFEOperator(tfeop_im_ad, tfeop_ex_ad)

tfeops = (
  tfeop_nl_man,
  tfeop_ql_man,
  tfeop_sl_man,
  tfeop_l_man,
  tfeop_imex_man,
  tfeop_nl_ad,
  tfeop_ql_ad,
  tfeop_sl_ad,
  tfeop_l_ad,
  tfeop_imex_ad,
)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
∂tuh0 = interpolate_everywhere(∂t(u)(t0), U0)

tol = 1.0e-6
sysslvr_l = LUSolver()
sysslvr_nl = NLSolver(sysslvr_l, show_trace=false, method=:newton, iterations=10)

# Testing function
function test_tfeop_order2(odeslvr, tfeop, uhs0)
  fesltn = solve(odeslvr, tfeop, t0, tF, uhs0)

  for (t_n, uh_n) in fesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
    @test e_n < tol
  end
end

odeslvrs = (
  Newmark(sysslvr_nl, dt, 0.5, 0.25),
)

uhs0 = (uh0, ∂tuh0)
for odeslvr in odeslvrs
  for tfeop in tfeops
    test_tfeop_order2(odeslvr, tfeop, uhs0)
  end
end

# Test with initial acceleration
∂ttuh0 = interpolate_everywhere(∂tt(u)(t0), U0)
uhs0 = (uh0, ∂tuh0, ∂ttuh0)
for odeslvr in odeslvrs
  for tfeop in tfeops
    test_tfeop_order2(odeslvr, tfeop, uhs0)
  end
end

end # module Order2FETests
