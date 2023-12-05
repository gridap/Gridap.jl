module Order2FETests

using Test

using LinearAlgebra
using ForwardDiff

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * (1 + t * (1 - t))
u(t::Real) = x -> u(x, t)
∂tu(t::Real) = ∂t(u)(t)
∂ttu(t::Real) = ∂tt(u)(t)

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
f(t) = x -> ∂tt(u)(x, t) + ∂t(u)(x, t) - Δ(u(t))(x)
mass(t, ∂ₜₜu, v) = ∫(∂ₜₜu ⋅ v) * dΩ
damping(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

res(t, u, v) = mass(t, ∂tt(u), v) + damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = damping(t, dut, v)
jac_tt(t, u, dutt, v) = mass(t, dutt, v)

mass_at_∂tt(t, u, v) = mass(t, ∂tt(u), v)
damping_at_∂t(t, u, v) = damping(t, ∂t(u), v)
res_quasilinear(t, u, v) = damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
res_linear(t, v) = (-1) * forcing(t, v)

feop_nonlinear = TransientFEOperator(res, jac, jac_t, jac_tt, U, V)
feop_quasilinear = TransientQuasilinearFEOperator(mass_at_∂tt, res_quasilinear, jac, jac_t, jac_tt, U, V)
feop_semilinear = TransientSemilinearFEOperator(mass_at_∂tt, res_quasilinear, jac, jac_t, jac_tt, U, V)
feop_linear = TransientLinearFEOperator(mass_at_∂tt, damping_at_∂t, stiffness, res_linear, jac, jac_t, jac_tt, U, V)

# Zero explicit residual
bilin0(t, u, v) = ∫(0 * u * v) * dΩ
jac0(t, u, du, v) = ∫(0 * du * v) * dΩ
feop_imex = TransientIMEXFEOperator(
  TransientSemilinearFEOperator(mass_at_∂tt, res_quasilinear, jac, jac_t, jac_tt, U, V),
  TransientFEOperator(bilin0, jac0, jac0, U, V)
)

feops = (
  feop_nonlinear,
  feop_quasilinear,
  feop_semilinear,
  feop_linear,
  feop_imex,
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

function test_transient_fe_problem_2(odeslvr, feop, uhs0)
  fesltn = solve(odeslvr, feop, uhs0, t0, tF)

  for (uh_n, t_n) in fesltn
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
    test_transient_fe_problem_2(odeslvr, feop, uhs0)
  end
end

# Solvers with memory
odeslvrs = (
  Newmark(disslvr, dt, 0.5, 0.25),
)

uhs0 = (uh0, ∂tuh0, ∂ttuh0,)
for odeslvr in odeslvrs
  for feop in feops
    test_transient_fe_problem_2(odeslvr, feop, uhs0)
  end
end

end # module Order2FETests
