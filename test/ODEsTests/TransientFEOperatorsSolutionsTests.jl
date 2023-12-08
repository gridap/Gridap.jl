module TransientFEOperatorsSolutionsTests

# This file tests both TransientFEOperators.jl and TransientFESolutions.jl

using Test

using LinearAlgebra
using LinearAlgebra: fillstored!

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs

include("ODESolversMocks.jl")

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

function test_transient_operator(
  feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms
)
  odeop = get_algebraic_operator(feop)
  odeopcache = allocate_odeopcache(odeop, t0, uhs0_dof)
  update_odeopcache!(odeopcache, odeop, t0)

  if !(feop isa TransientIMEXFEOperator)
    @test test_transient_fe_operator(feop, t0, uhs0_cf)
  end

  fesltn = solve(odeslvr, feop, t0, tF, uhs0_fe)
  @test test_transient_fe_solution(fesltn)

  # Test storage of constant forms
  if feop isa TransientIMEXFEOperator
    im_feop, ex_feop = feop.im_feop, feop.ex_feop
    im_odeopcache, ex_odeopcache, _ = odeopcache
    feops_odeopcaches = ((im_feop, im_odeopcache), (ex_feop, ex_odeopcache))
  else
    feops_odeopcaches = ((feop, odeopcache),)
  end
  for (feop, odeopcache) in feops_odeopcaches
    forms = get_forms(feop)
    for k in 0:length(forms)-1
      constant_form = constant_forms[k+1]
      @test is_form_constant(feop, k) == constant_form
      if constant_form
        @test !isnothing(odeopcache.const_forms[k+1])
      end
    end
  end
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

# Initial data
uhs0_fe = (uh0,)
uhs0_cf = TransientCellField(uh0, (∂ₜuh0,))
uhs0_dof = (uh0_dof, ∂ₜuh0_dof)

# TransientFEOperator
constant_forms = ()

feop = TransientFEOperator(res, (jac, jac_t), U, V)
test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

feop = TransientFEOperator(res, U, V; order)
test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

# TransientQuasilinearFEOperator
constant_forms = (false,)

feop = TransientQuasilinearFEOperator(mass, res_ql, (jac, jac_t), U, V)
test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

feop = TransientQuasilinearFEOperator(mass, res_ql, U, V; order)
test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

for constant_mass in (true, false)
  # TransientSemilinearFEOperator
  constant_forms = (constant_mass,)

  feop = TransientSemilinearFEOperator(
    mass, res_ql, (jac, jac_t), U, V;
    constant_mass
  )
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  feop = TransientSemilinearFEOperator(
    mass, res_ql, U, V;
    constant_mass, order
  )
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  # TransientIMEXFEOperator
  constant_forms = (constant_mass,)

  im_feop = TransientSemilinearFEOperator(
    mass, im_res, (im_jac, im_jac_t), U, V;
    constant_mass
  )
  ex_feop = TransientFEOperator(ex_res, (ex_jac,), U, V)
  feop = TransientIMEXFEOperator(im_feop, ex_feop)
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  for constant_stiffness in (true, false)
    # TransientLinearFEOperator
    constant_forms = (constant_mass, constant_stiffness)

    feop = TransientLinearFEOperator(
      (stiffness, mass), res_l, (jac, jac_t), U, V;
      constant_forms
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

    feop = TransientLinearFEOperator(
      (stiffness, mass), res_l, U, V;
      constant_forms
    )
    test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)
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

# Initial data
uhs0_fe = (uh0, ∂ₜuh0)
uhs0_cf = TransientCellField(uh0, (∂ₜuh0, ∂ₜₜuh0,))
uhs0_dof = (uh0_dof, ∂ₜuh0_dof, ∂ₜₜuh0_dof)

# TransientFEOperator
constant_forms = ()

feop = TransientFEOperator(res, (jac, jac_t, jac_tt), U, V)
test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

feop = TransientFEOperator(res, U, V; order)
test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

for constant_mass in (true, false)
  # TransientQuasilinearFEOperator
  constant_forms = (false,)

  feop = TransientQuasilinearFEOperator(mass, res_ql, (jac, jac_t, jac_tt), U, V)
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  feop = TransientQuasilinearFEOperator(mass, res_ql, U, V; order)
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  # TransientSemilinearFEOperator
  constant_forms = (constant_mass,)

  feop = TransientSemilinearFEOperator(
    mass, res_ql, (jac, jac_t, jac_tt), U, V;
    constant_mass
  )
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  feop = TransientSemilinearFEOperator(
    mass, res_ql, U, V;
    constant_mass, order
  )
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  # TransientIMEXFEOperator
  constant_forms = (constant_mass,)

  im_feop = TransientSemilinearFEOperator(
    mass, im_res, (im_jac, im_jac_t, im_jac_tt), U, V;
    constant_mass
  )
  ex_feop = TransientFEOperator(ex_res, (ex_jac, ex_jac_t), U, V)
  feop = TransientIMEXFEOperator(im_feop, ex_feop)
  test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  for constant_damping in (true, false)
    for constant_stiffness in (true, false)
      # TransientLinearFEOperator
      constant_forms = (constant_stiffness, constant_damping, constant_mass)

      feop = TransientLinearFEOperator(
        (stiffness, damping, mass), res_l, (jac, jac_t, jac_tt), U, V;
        constant_forms
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

      feop = TransientLinearFEOperator(
        (stiffness, damping, mass), res_l, U, V;
        constant_forms
      )
      test_transient_operator(feop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)
    end
  end
end

end # module TransientFEOperatorsSolutionsTests
