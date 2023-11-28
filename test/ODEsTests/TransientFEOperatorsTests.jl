module TransientFEOperatorsTests

using Test

using LinearAlgebra
using LinearAlgebra: fillstored!
using ForwardDiff

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

# Analytical functions
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * (t + 3.0)
u(t::Real) = x -> u(x, t)
u(x) = t -> u(x, t)

∂tu(x, t) = ∂t(u)(x, t)
∂tu(t::Real) = x -> ∂tu(x, t)

f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)

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

# ODE operator
m(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
a(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
b(t, v) = ∫(f(t) ⋅ v) * dΩ

t0 = 0.0
t1 = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)

# Residual and jacobian with FEOperator
dt⁻¹ = inv(dt)
_res(u, v) = dt⁻¹ * m(t0, u, v) + a(t0, u, v) - b(t0, v)
_jac(u, du, v) = dt⁻¹ * m(t0, du, v) + a(t0, du, v)
_op = FEOperator(_res, _jac, U0, V)

_r = residual(_op, uh0)
_J = jacobian(_op, uh0)

# Residual and jacobian with TransientFEOperators
# Testing with all combinations of constant jacobians and forcing term,
# With manual or automatic jacobians
∂ₜuh0 = FEFunction(U0, get_free_dof_values(uh0) ./ dt)
uₜ = TransientCellField(uh0, (∂ₜuh0,))
us = (get_free_dof_values(uh0), get_free_dof_values(∂ₜuh0))

mass(t, u, v) = m(t, ∂t(u), v)
stiffness(t, u, v) = a(t, u, v)
res(t, u, v) = mass(t, u, v) + stiffness(t, u, v) - b(t, v)
jac(t, u, du, v) = a(t, du, v)
jac_t(t, u, dut, v) = m(t, dut, v)

res_masslinear(t, u, v) = a(t, u, v) - b(t, v)
res_linear(t, v) = (-1) * b(t, v)

function test_transient_operator(op)
  @test test_transient_fe_operator(op, t0, uₜ)

  ode_op = get_algebraic_operator(op)
  ode_cache = allocate_cache(ode_op, t0, us)

  r = allocate_residual(ode_op, t0, us, ode_cache)
  J = allocate_jacobian(ode_op, t0, us, ode_cache)
  residual!(r, ode_op, t0, us, ode_cache)
  @test all(r .≈ _r)

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, ode_op, t1, us, (1, dt⁻¹), ode_cache)
  @test all(J .≈ _J)
end

for jac_u_constant in (true, false)
  for jac_u̇_constant in (true, false)
    jacs_constant = (jac_u_constant, jac_u̇_constant)

    # TransientFEOperator
    op = TransientFEOperator(res, jac, jac_t, U, V; jacs_constant)
    test_transient_operator(op)

    op = TransientFEOperator(res, U, V; jacs_constant)
    test_transient_operator(op)

    # TransientMassLinearFEOperator
    op = TransientMassLinearFEOperator(
      mass, res_masslinear, jac, jac_t, U, V;
      jacs_constant, forcing_constant=false
    )
    test_transient_operator(op)

    op = TransientMassLinearFEOperator(
      mass, res_masslinear, U, V;
      jacs_constant, forcing_constant=false
    )
    test_transient_operator(op)

    # TransientLinearFEOperator
    op = TransientLinearFEOperator(
      mass, stiffness, res_linear, jac, jac_t, U, V;
      jacs_constant, forcing_constant=false
    )
    test_transient_operator(op)

    op = TransientLinearFEOperator(
      (mass, stiffness), res_linear, U, V;
      jacs_constant, forcing_constant=false
    )
    test_transient_operator(op)
  end
end

# TODO Test second-order operators

end # module TransientFEOperatorsTests
