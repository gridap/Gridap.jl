module FreeSurfacePotentialFlowTests

using Test

using Gridap
using Gridap.Geometry

# Parameters
L = 2 * π
H = 1.0
n = 8
order = 2
g = 9.81
ξ = 0.1
λ = L / 2
k = 2 * π / L
h = L / n
ω = √(g * k * tanh(k * H))
t0 = 0.0
dt = h / (2 * λ * ω)
tF = 10 * dt # 2 * π
α = 2 / dt
tol = 1.0e-2

# Exact solution
ϕₑ(t) = x -> ω / k * ξ * (cosh(k * (x[2]))) / sinh(k * H) * sin(k * x[1] - ω * t)
ηₑ(t) = x -> ξ * cos(k * x[1] - ω * t)

# Domain
domain = (0, L, 0, H)
partition = (n, n)
model = CartesianDiscreteModel(domain, partition; isperiodic=(true, false))

# Boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "bottom", [1, 2, 5])
add_tag_from_tags!(labels, "free_surface", [3, 4, 6])

# Triangulation
Ω = Interior(model)
Γ = Boundary(model, tags="free_surface")
dΩ = Measure(Ω, 2 * order)
dΓ = Measure(Γ, 2 * order)

# FE spaces
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(Ω, reffe, conformity=:H1)
V_Γ = TestFESpace(Γ, reffe, conformity=:H1)
U = TransientTrialFESpace(V)
U_Γ = TransientTrialFESpace(V_Γ)
X = TransientMultiFieldFESpace([U, U_Γ])
Y = MultiFieldFESpace([V, V_Γ])

# Weak form
m(t, (ϕt, ηt), (w, v)) = ∫(0.5 * (α / g * (w * ϕt) + v * ϕt) - (w * ηt))dΓ
a(t, (ϕ, η), (w, v)) = ∫(∇(ϕ) ⋅ ∇(w))dΩ + ∫(0.5 * (α * (w * η) + g * v * η))dΓ
b(t, (w, v)) = ∫(0.0 * w)dΓ

res(t, x, y) = m(t, ∂t(x), y) + a(t, x, y) - b(t, y)
jac(t, x, dx, y) = a(t, dx, y)
jac_t(t, x, dxt, y) = m(t, dxt, y)

# Optimal transient FE Operator
op_const = TransientLinearFEOperator((a, m), b, (jac, jac_t), X, Y, constant_forms=(true, true))

# TransientFEOperator exploiting automatic differentiation (testing purposes)
op_trans = TransientFEOperator(res, (jac, jac_t), X, Y)
op_ad = TransientFEOperator(res, X, Y)

# TransientFEOperator exploiting time derivative of separate fields (TransientMultiFieldCellField)
res2(t, (ϕ, η), y) = m(t, (∂t(ϕ), ∂t(η)), y) + a(t, (ϕ, η), y) - b(t, y)
op_multifield = TransientFEOperator(res2, (jac, jac_t), X, Y)

# Solver
sysslvr_l = LUSolver()
sysslvr_nl = NLSolver(sysslvr_l, show_trace=false, method=:newton, iterations=10)
odeslvr = ThetaMethod(sysslvr_nl, dt, 0.5)

# Initial solution
U0 = U(t0)
UΓ0 = U_Γ(t0)
X0 = X(t0)
uh0 = interpolate_everywhere(ϕₑ(t0), U0)
uhΓ0 = interpolate_everywhere(ηₑ(t0), UΓ0)
xh0 = interpolate_everywhere([uh0, uhΓ0], X0)
xhs0 = (xh0,)

function test_flow_operator(op)
  fesltn = solve(odeslvr, op, t0, tF, xhs0)

  # Post-process
  l2_Ω(v) = √(∑(∫(v ⋅ v) * dΩ))
  l2_Γ(v) = √(∑(∫(v ⋅ v) * dΓ))
  E_kin(v) = 0.5 * ∑(∫(∇(v) ⋅ ∇(v)) * dΩ)
  E_pot(v) = g * 0.5 * ∑(∫(v * v)dΓ)
  Eₑ = 0.5 * g * ξ^2 * L

  for (tn, (ϕn, ηn)) in fesltn
    E = E_kin(ϕn) + E_pot(ηn)
    error_ϕ = l2_Ω(ϕn - ϕₑ(tn))
    error_η = l2_Γ(ηn - ηₑ(tn))
    @test abs(E / Eₑ - 1.0) <= tol
    @test error_ϕ <= tol
    @test error_η <= tol
  end
end

# op_ad not working yet
for op in (op_const, op_trans, op_multifield)
  test_flow_operator(op)
end

end
