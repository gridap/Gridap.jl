module HeatEquationNeumannTests

using Test

using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
u(x, t) = x[1] * (1 - x[2]) * (1 + t)
u(t::Real) = x -> u(x, t)
u(x) = t -> u(x, t)

∂tu(x, t) = ∂t(u)(x, t)
∂tu(t::Real) = x -> ∂tu(x, t)

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
f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)
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
res_l(t, v) = (-1) * forcing(t, v)

args_man = ((jac, jac_t), U, V)
feop_nl_man = TransientFEOperator(res, args_man...)
feop_ql_man = TransientQuasilinearFEOperator(mass, res_ql, args_man...)
feop_sl_man = TransientSemilinearFEOperator(mass, res_ql, args_man...)
feop_l_man = TransientLinearFEOperator((stiffness, mass), res_l, args_man...)

args_ad = (U, V)
feop_nl_ad = TransientFEOperator(res, args_ad...)
feop_ql_ad = TransientQuasilinearFEOperator(mass, res_ql, args_ad...)
feop_sl_ad = TransientSemilinearFEOperator(mass, res_ql, args_ad...)
feop_l_ad = TransientLinearFEOperator((stiffness, mass), res_l, args_ad...)

feops = (
  feop_nl_man,
  feop_ql_man,
  feop_sl_man,
  feop_l_man,
  feop_nl_ad,
  feop_ql_ad,
  feop_sl_ad,
  feop_l_ad
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
disslvr = LUSolver()
odeslvrs = (
  MidPoint(disslvr, dt),
)

# Tests
for odeslvr in odeslvrs
  for feop in feops
    fesltn = solve(odeslvr, feop, t0, tF, uhs0)

    for (t_n, uh_n) in fesltn
      eh_n = u(t_n) - uh_n
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol
    end
  end
end

end # module HeatEquationNeumannTests
