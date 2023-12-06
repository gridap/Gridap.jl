module TransientFEOperatorsSolutionsTests

# This files tests both TransientFEOperators.jl and TransientFESolutions.jl

using Test

using LinearAlgebra
using LinearAlgebra: fillstored!

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

include("ODESolversMocks.jl")

# Analytical functions
u(x, t) = (1.0 - x[1]) * x[1] * (1.0 - x[2]) * x[2] * (1 + t * (1 - t))
u(t::Real) = x -> u(x, t)
u(x) = t -> u(x, t)

∂tu(x, t) = ∂t(u)(x, t)
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

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.2
dt⁻¹ = 1 / dt
dt⁻² = dt⁻¹^2

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
∂ₜuh0 = FEFunction(U0, get_free_dof_values(uh0) .* dt⁻¹)
∂ₜₜuh0 = FEFunction(U0, get_free_dof_values(uh0) .* dt⁻²)

uh0_dof = get_free_dof_values(uh0)
∂ₜuh0_dof = get_free_dof_values(∂ₜuh0)
∂ₜₜuh0_dof = get_free_dof_values(∂ₜₜuh0)

# ODE solver
atol = 1.0e-12
rtol = 1.0e-8
maxiter = 100
disslvr = DiscreteODESolverMock(rtol, atol, maxiter)
odeslvr = ODESolverMock(disslvr, dt)

function test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)
  odeop = get_algebraic_operator(feop)
  odeopcache = allocate_odeopcache(odeop, t0, uhs0_dof)

  r = allocate_residual(odeop, t0, uhs0_dof, odeopcache)
  J = allocate_jacobian(odeop, t0, uhs0_dof, odeopcache)
  residual!(r, odeop, t0, uhs0_dof, odeopcache)
  @test r ≈ r_const

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, odeop, t0, uhs0_dof, ws, odeopcache)
  @test J ≈ J_const

  if !(feop isa TransientIMEXFEOperator)
    @test test_transient_fe_operator(feop, t0, uhs0_cf)
  end

  fesltn = solve(odeslvr, feop, uhs0_fe, t0, tF)
  @test test_transient_fe_solution(fesltn)
end

###############
# First-order #
###############
order = 1
f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)

mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
mass(t, u, ∂ₜu, v) = mass(t, ∂ₜu, v)
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

# Residual and jacobian with FEOperator
# The TransientFEOperator will be evaluated at (uh0, uh0 / dt) so we add
# the factor dt⁻¹ in the FEOperator for comparison
res_const(u, v) = dt⁻¹ * mass(t0, u, v) + stiffness(t0, u, v) - forcing(t0, v)
jac_const(u, du, v) = dt⁻¹ * mass(t0, u, du, v) + stiffness(t0, du, v)
feop_const = FEOperator(res_const, jac_const, U0, V)
r_const = residual(feop_const, uh0)
J_const = jacobian(feop_const, uh0)

# Residual and jacobian with TransientFEOperators
# Testing with all combinations of constant jacobians and residual,
# With manual or automatic jacobians
uhs0_fe = (uh0,)
uhs0_cf = TransientCellField(uh0, (∂ₜuh0,))
uhs0_dof = (uh0_dof, ∂ₜuh0_dof)
ws = (1, dt⁻¹)

res(t, u, v) = mass(t, u, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = mass(t, u, dut, v)

res_ql(t, u, v) = stiffness(t, u, v) - forcing(t, v)
res_l(t, v) = (-1) * forcing(t, v)

# TODO could think of a simple and optimised way to create a zero residual or
# jacobian without assembling the vector / matrix
im_res(t, u, v) = ∫(0 * u * v) * dΩ
im_jac(t, u, du, v) = ∫(0 * du * v) * dΩ
im_jac_t(t, u, dut, v) = mass(t, u, dut, v)

ex_res(t, u, v) = stiffness(t, u, v) - forcing(t, v)
ex_jac(t, u, du, v) = stiffness(t, du, v)

for jac_u_constant in (true, false)
  for jac_u̇_constant in (true, false)
    jacs_constant = (jac_u_constant, jac_u̇_constant)
    residual_constant = false

    # TransientFEOperator
    feop = TransientFEOperator(
      res, (jac, jac_t), U, V;
      jacs_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    feop = TransientFEOperator(
      res, U, V;
      jacs_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    # TransientQuasilinearFEOperator
    feop = TransientQuasilinearFEOperator(
      mass, res_ql, (jac, jac_t), U, V;
      jacs_constant, residual_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    feop = TransientQuasilinearFEOperator(
      mass, res_ql, U, V;
      jacs_constant, residual_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    # TransientSemilinearFEOperator
    feop = TransientSemilinearFEOperator(
      mass, res_ql, (jac, jac_t), U, V;
      jacs_constant, residual_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    feop = TransientSemilinearFEOperator(
      mass, res_ql, U, V;
      jacs_constant, residual_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    # TransientLinearFEOperator
    feop = TransientLinearFEOperator(
      (stiffness, mass), res_l, (jac, jac_t), U, V;
      jacs_constant, residual_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    feop = TransientLinearFEOperator(
      (stiffness, mass), res_l, U, V;
      jacs_constant, residual_constant
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

    # TransientIMEXFEOperator
    im_feop = TransientSemilinearFEOperator(
      mass, im_res, (im_jac, im_jac_t), U, V;
      jacs_constant, residual_constant
    )
    ex_feop = TransientFEOperator(
      ex_res, (ex_jac,), U, V;
      jacs_constant=(jac_u_constant,)
    )
    feop = TransientIMEXFEOperator(im_feop, ex_feop)
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)
  end
end

################
# Second-order #
################
order = 2
f(t) = x -> ∂tt(u)(x, t) + ∂t(u)(x, t) - Δ(u(t))(x)

mass(t, ∂ₜₜu, v) = ∫(∂ₜₜu ⋅ v) * dΩ
mass(t, u, ∂ₜₜu, v) = mass(t, ∂ₜₜu, v)
damping(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

# Residual and jacobian with FEOperator
dt⁻² = dt⁻¹^2
res_const(u, v) = dt⁻² * mass(t0, u, v) + dt⁻¹ * damping(t0, u, v) + stiffness(t0, u, v) - forcing(t0, v)
jac_const(u, du, v) = dt⁻² * mass(t0, du, v) + dt⁻¹ * damping(t0, du, v) + stiffness(t0, du, v)
feop_const = FEOperator(res_const, jac_const, U0, V)
r_const = residual(feop_const, uh0)
J_const = jacobian(feop_const, uh0)

# Residual and jacobian with TransientFEOperators
uhs0_fe = (uh0, ∂ₜuh0)
uhs0_cf = TransientCellField(uh0, (∂ₜuh0, ∂ₜₜuh0,))
uhs0_dof = (uh0_dof, ∂ₜuh0_dof, ∂ₜₜuh0_dof)
ws = (1, dt⁻¹, dt⁻²)

res(t, u, v) = mass(t, u, ∂tt(u), v) + damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = damping(t, dut, v)
jac_tt(t, u, dutt, v) = mass(t, dutt, v)

res_ql(t, u, v) = damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
res_l(t, v) = (-1) * forcing(t, v)

im_res(t, u, v) = ∫(0 * u * v) * dΩ
im_jac(t, u, du, v) = ∫(0 * du * v) * dΩ
im_jac_t(t, u, dut, v) = ∫(0 * dut * v) * dΩ
im_jac_tt(t, u, dutt, v) = mass(t, u, dutt, v)

ex_res(t, u, v) = damping(t, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
ex_jac(t, u, du, v) = stiffness(t, du, v)
ex_jac_t(t, u, dut, v) = damping(t, dut, v)

for jac_u_constant in (true, false)
  for jac_u̇_constant in (true, false)
    for jac_ü_constant in (true, false)
      jacs_constant = (jac_u_constant, jac_u̇_constant, jac_ü_constant)
      residual_constant = false

      # TransientFEOperator
      feop = TransientFEOperator(
        res, (jac, jac_t, jac_tt), U, V;
        jacs_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      feop = TransientFEOperator(
        res, U, V;
        order, jacs_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      # TransientQuasilinearFEOperator
      feop = TransientQuasilinearFEOperator(
        mass, res_ql, (jac, jac_t, jac_tt), U, V;
        jacs_constant, residual_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      feop = TransientQuasilinearFEOperator(
        mass, res_ql, U, V;
        order, jacs_constant, residual_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      # TransientSemilinearFEOperator
      feop = TransientSemilinearFEOperator(
        mass, res_ql, (jac, jac_t, jac_tt), U, V;
        jacs_constant, residual_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      feop = TransientSemilinearFEOperator(
        mass, res_ql, U, V;
        order, jacs_constant, residual_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      # TransientLinearFEOperator
      feop = TransientLinearFEOperator(
        (stiffness, damping, mass), res_l, (jac, jac_t, jac_tt), U, V;
        jacs_constant, residual_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      feop = TransientLinearFEOperator(
        (stiffness, damping, mass), res_l, U, V;
        order, jacs_constant, residual_constant
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)

      # TransientIMEXFEOperator
      im_feop = TransientSemilinearFEOperator(
        mass, im_res, (im_jac, im_jac_t, im_jac_tt), U, V;
        jacs_constant, residual_constant
      )
      ex_feop = TransientFEOperator(
        ex_res, (ex_jac, ex_jac_t), U, V;
        jacs_constant=(jac_u_constant, jac_u̇_constant)
      )
      feop = TransientIMEXFEOperator(im_feop, ex_feop)
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, ws)
    end
  end
end

end # module TransientFEOperatorsSolutionsTests
