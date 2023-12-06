module HeatEquationNeumannTests

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
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing_Ω(t, v) = ∫(f(t) ⋅ v) * dΩ
forcing_Γ(t, v) = ∫((∇(u(t)) ⋅ nΓ) ⋅ v) * dΓ
forcing(t, v) = forcing_Ω(t, v) + forcing_Γ(t, v)

res(t, u, v) = mass(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = mass(t, dut, v)

feop = TransientFEOperator(res, jac, jac_t, U, V)
feop_AD = TransientFEOperator(res, U, V)
feops = (feop, feop_AD,)

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
    fesltn = solve(odeslvr, feop, uhs0, t0, tF)

    for (uh_n, t_n) in fesltn
      eh_n = u(t_n) - uh_n
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol
    end
  end
end

end # module HeatEquationNeumannTests
