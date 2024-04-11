module HeatEquationDGTests

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

# FE spaces
order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V = FESpace(model, reffe, conformity=:L2)
U = TransientTrialFESpace(V, u)

# Integration
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ, degree)
nΓ = get_normal_vector(Γ)

Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ, degree)
nΛ = get_normal_vector(Λ)

# FE operator
ft(t) = x -> ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)
h = 1 / 5
γ = order * (order + 1)

mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
stiffness_Ω(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
stiffness_Γ(t, u, v) = ∫((γ / h) * (u ⋅ v) - (∇(u) ⋅ nΓ) ⋅ v - u ⋅ (∇(v) ⋅ nΓ)) * dΓ
stiffness_Λ(t, u, v) = ∫((γ / h) * (jump(u * nΛ) ⊙ jump(v * nΛ)) - mean(∇(u)) ⊙ jump(v * nΛ) - jump(u * nΛ) ⊙ mean(∇(v))) * dΛ
stiffness(t, u, v) = stiffness_Ω(t, u, v) + stiffness_Γ(t, u, v) + stiffness_Λ(t, u, v)
forcing_Ω(t, v) = ∫(f(t) ⋅ v) * dΩ
forcing_Γ(t, v) = ∫((γ / h) * v ⋅ u(t) - u(t) ⋅ (∇(v) ⋅ nΓ)) * dΓ
forcing(t, v) = forcing_Ω(t, v) + forcing_Γ(t, v)

res(t, u, v) = mass(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = mass(t, dut, v)

tfeop_nl_man = TransientFEOperator(res, (jac, jac_t), U, V)

# TODO there is an issue with AD here. The issue is already there in the
# current version of Gridap.ODEs. This happens when calling
# TransientCellFieldType(y, u.derivatives) in the construction of the jacobians
# with AD
tfeop_nl_ad = TransientFEOperator(res, U, V)

tfeops = (
  tfeop_nl_man,
  # tfeop_nl_ad
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

end # module HeatEquationDGTests
