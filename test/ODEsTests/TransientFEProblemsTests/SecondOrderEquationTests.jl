module Order2FETests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
u(x, t) = x[1] * (1 - x[2]) * (1 + t)
∂tu(x, t) = ∂t(u)(x, t)
∂ttu(x, t) = ∂tt(u)(x, t)

u(t::Real) = x -> u(x, t)
∂tu(t::Real) = x -> ∂tu(x, t)
∂ttu(t::Real) = x -> ∂ttu(x, t)

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
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# FE operator
order = 2
f(t) = x -> ∂tt(u)(x, t) + ∂t(u)(x, t) - Δ(u(t))(x)

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
res_l(t, v) = (-1) * forcing(t, v)

res0(t, u, v) = ∫(0 * u * v) * dΩ
jac0(t, u, du, v) = ∫(0 * du * v) * dΩ

args_man = ((jac, jac_t, jac_tt), U, V)
args0_man = ((jac0, jac0), U, V)
feop_nl_man = TransientFEOperator(res, args_man...)
feop_ql_man = TransientQuasilinearFEOperator(mass, res_ql, args_man...)
feop_sl_man = TransientSemilinearFEOperator(mass, res_ql, args_man...)
feop_l_man = TransientLinearFEOperator((stiffness, damping, mass), res_l, args_man...)

feop_im_man = TransientSemilinearFEOperator(mass, res_ql, args_man...)
feop_ex_man = TransientFEOperator(res0, args0_man...)
feop_imex_man = TransientIMEXFEOperator(feop_im_man, feop_ex_man)

args_ad = (U, V)
feop_nl_ad = TransientFEOperator(res, args_ad...)
feop_ql_ad = TransientQuasilinearFEOperator(mass, res_ql, args_ad...)
feop_sl_ad = TransientSemilinearFEOperator(mass, res_ql, args_ad...)
feop_l_ad = TransientLinearFEOperator((stiffness, damping, mass), res_l, args_ad...)

feop_im_ad = TransientSemilinearFEOperator(mass, res_ql, args_ad...)
feop_ex_ad = TransientFEOperator(res0, args_ad..., order=0)
feop_imex_ad = TransientIMEXFEOperator(feop_im_ad, feop_ex_ad)

feops = (
  feop_nl_man,
  feop_ql_man,
  feop_sl_man,
  feop_l_man,
  feop_imex_man,
  # feop_nl_ad,
  # feop_ql_ad,
  # feop_sl_ad,
  feop_l_ad,
  # feop_imex_ad,
)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
∂tuh0 = interpolate_everywhere(∂tu(t0), U0)
∂ttuh0 = interpolate_everywhere(∂ttu(t0), U0)

# Allow for slightly larger errors because the forcing term is second-order in
# time here (so that the second-derivative is not zero), compared to
# first-order in Order1Tests.jl
tol = 1.0e-5
disslvr = LUSolver()

# Testing function
function test_feop_order2(odeslvr, feop, uhs0)
  fesltn = solve(odeslvr, feop, t0, tF, uhs0)

  for (t_n, uh_n) in fesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
    @test e_n < tol
  end
end

# Solvers without memory
odeslvrs = ()

uhs0 = (uh0, ∂tuh0,)
for odeslvr in odeslvrs
  for feop in feops
    test_feop_order2(odeslvr, feop, uhs0)
  end
end

# Solvers with memory
odeslvrs = (
  Newmark(disslvr, dt, 0.5, 0.25),
)

uhs0 = (uh0, ∂tuh0, ∂ttuh0,)
for odeslvr in odeslvrs
  for feop in feops
    test_feop_order2(odeslvr, feop, uhs0)
  end
end

end # module Order2FETests
