module Order1FETests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * (1 + t)
∂tu(x, t) = ∂t(u)(x, t)

u(t::Real) = x -> u(x, t)
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

# Integration
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# FE operator
f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)

mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
mass(t, u, ∂ₜu, v) = mass(t, ∂ₜu, v)
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

res(t, u, v) = mass(t, u, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = mass(t, u, dut, v)

res_ql(t, u, v) = stiffness(t, u, v) - forcing(t, v)
res_l(t, v) = (-1) * forcing(t, v)

res0(t, u, v) = ∫(0 * u * v) * dΩ
jac0(t, u, du, v) = ∫(0 * du * v) * dΩ

args = ((jac, jac_t), U, V)
args0 = ((jac0,), U, V)
feop_nl = TransientFEOperator(res, args...)
feop_ql = TransientQuasilinearFEOperator(mass, res_ql, args...)
feop_sl = TransientSemilinearFEOperator(mass, res_ql, args...)
feop_l = TransientLinearFEOperator((stiffness, mass), res_l, args...)

feop_im = TransientSemilinearFEOperator(mass, res_ql, args...)
feop_ex = TransientFEOperator(res0, args0...)
feop_imex = TransientIMEXFEOperator(feop_im, feop_ex)

feops = (
  feop_nl,
  feop_ql,
  feop_sl,
  feop_l,
  feop_imex,
)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
∂tuh0 = interpolate_everywhere(∂tu(t0), U0)

tol = 1.0e-6
disslvr = LUSolver()

# Testing function
function test_feop_order1(odeslvr, feop, uhs0)
  fesltn = solve(odeslvr, feop, uhs0, t0, tF)

  for (uh_n, t_n) in fesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
    @test e_n < tol
  end
end

# Solvers without memory
# Do not try explicit solvers
odeslvrs = (
  ThetaMethod(disslvr, dt, 0.2),
  MidPoint(disslvr, dt),
  ThetaMethod(disslvr, dt, 0.8),
  BackwardEuler(disslvr, dt),
  RungeKutta(disslvr, disslvr, dt, :BE_1_0_1),
  RungeKutta(disslvr, disslvr, dt, :CN_2_0_2),
  RungeKutta(disslvr, disslvr, dt, :SDIRK_2_0_2),
  RungeKutta(disslvr, disslvr, dt, :SDIRK_2_0_3),
  RungeKutta(disslvr, disslvr, dt, :ESDIRK_3_1_2),
  RungeKutta(disslvr, disslvr, dt, :TRBDF2_3_2_3),
  RungeKutta(disslvr, disslvr, dt, :TRX2_3_2_3),
)

uhs0 = (uh0,)
for odeslvr in odeslvrs
  for feop in feops
    test_feop_order1(odeslvr, feop, uhs0)
  end
end

# Solvers with memory
odeslvrs = (
  GeneralizedAlpha1(disslvr, dt, 0.0),
  GeneralizedAlpha1(disslvr, dt, 0.5),
  GeneralizedAlpha1(disslvr, dt, 1.0),
)

uhs0 = (uh0, ∂tuh0)
for odeslvr in odeslvrs
  for feop in feops
    test_feop_order1(odeslvr, feop, uhs0)
  end
end

# Solvers for IMEX decompositions
feops = (
  feop_imex,
)

odeslvrs = (
  RungeKutta(disslvr, disslvr, dt, :IMEX_Midpoint_2_0_2),
)

uhs0 = (uh0,)
for odeslvr in odeslvrs
  for feop in feops
    test_feop_order1(odeslvr, feop, uhs0)
  end
end

end # module Order1FETests
