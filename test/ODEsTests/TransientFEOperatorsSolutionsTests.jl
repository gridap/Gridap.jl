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
ut(t) = x -> x[1] * (1 - x[2]) * (1 + t)
u = TimeSpaceFunction(ut)

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
sysslvr = NonlinearSolverMock(rtol, atol, maxiter)
odeslvr = ODESolverMock(sysslvr, dt)

function test_transient_operator(
  tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms
)
  odeop = get_algebraic_operator(tfeop)
  odeopcache = allocate_odeopcache(odeop, t0, uhs0_dof)
  update_odeopcache!(odeopcache, odeop, t0)

  if !(tfeop isa TransientIMEXFEOperator)
    @test test_tfe_operator(tfeop, t0, uhs0_cf)
  end

  fesltn = solve(odeslvr, tfeop, t0, tF, uhs0_fe)
  @test test_tfe_solution(fesltn)

  # Test storage of constant forms
  if tfeop isa TransientIMEXFEOperator
    im_tfeop, ex_tfeop = get_imex_operators(tfeop)
    im_odeopcache, ex_odeopcache = odeopcache
    tfeops_odeopcaches = ((im_tfeop, im_odeopcache), (ex_tfeop, ex_odeopcache))
  else
    tfeops_odeopcaches = ((tfeop, odeopcache),)
  end
  for (tfeop, odeopcache) in tfeops_odeopcaches
    num_forms = get_num_forms(tfeop)
    if num_forms == 1
      order = get_order(tfeop)
      constant_form = constant_forms[1]
      @test is_form_constant(tfeop, order) == constant_form
      if constant_form
        @test !isnothing(odeopcache.const_forms[1])
      end
    else
      for k in 0:num_forms-1
        constant_form = constant_forms[k+1]
        @test is_form_constant(tfeop, k) == constant_form
        if constant_form
          @test !isnothing(odeopcache.const_forms[k+1])
        end
      end
    end
  end
end

###############
# First-order #
###############
order = 1
ft(t) = x -> ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)

mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
mass(t, u, ∂ₜu, v) = mass(t, ∂ₜu, v)
stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

res(t, u, v) = mass(t, u, ∂t(u), v) + stiffness(t, u, v) - forcing(t, v)
jac(t, u, du, v) = stiffness(t, du, v)
jac_t(t, u, dut, v) = mass(t, u, dut, v)

res_ql(t, u, v) = stiffness(t, u, v) - forcing(t, v)
res_l(t, v) = forcing(t, v)

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

tfeop = TransientFEOperator(res, (jac, jac_t), U, V)
test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

tfeop = TransientFEOperator(res, U, V; order)
test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

# TransientQuasilinearFEOperator
constant_forms = (false,)

tfeop = TransientQuasilinearFEOperator(mass, res_ql, (jac, jac_t), U, V)
test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

tfeop = TransientQuasilinearFEOperator(mass, res_ql, U, V; order)
test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

for constant_mass in (true, false)
  # TransientSemilinearFEOperator
  constant_forms = (constant_mass,)

  tfeop = TransientSemilinearFEOperator(
    mass, res_ql, (jac, jac_t), U, V;
    constant_mass
  )
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  tfeop = TransientSemilinearFEOperator(
    mass, res_ql, U, V;
    constant_mass, order
  )
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  # TransientIMEXFEOperator
  constant_forms = (constant_mass,)

  im_tfeop = TransientSemilinearFEOperator(
    mass, im_res, (im_jac, im_jac_t), U, V;
    constant_mass
  )
  ex_tfeop = TransientFEOperator(ex_res, (ex_jac,), U, V)
  tfeop = TransientIMEXFEOperator(im_tfeop, ex_tfeop)
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  for constant_stiffness in (true, false)
    # TransientLinearFEOperator
    constant_forms = (constant_stiffness, constant_mass)

    tfeop = TransientLinearFEOperator(
      (stiffness, mass), res_l, (jac, jac_t), U, V;
      constant_forms
    )
    test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

    tfeop = TransientLinearFEOperator(
      (stiffness, mass), res_l, U, V;
      constant_forms
    )
    test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)
  end
end

################
# Second-order #
################
order = 2
ft(t) = x -> ∂tt(u)(t, x) + ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)

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
res_l(t, v) = forcing(t, v)

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

tfeop = TransientFEOperator(res, (jac, jac_t, jac_tt), U, V)
test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

tfeop = TransientFEOperator(res, U, V; order)
test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

for constant_mass in (true, false)
  # TransientQuasilinearFEOperator
  constant_forms = (false,)

  tfeop = TransientQuasilinearFEOperator(mass, res_ql, (jac, jac_t, jac_tt), U, V)
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  tfeop = TransientQuasilinearFEOperator(mass, res_ql, U, V; order)
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  # TransientSemilinearFEOperator
  constant_forms = (constant_mass,)

  tfeop = TransientSemilinearFEOperator(
    mass, res_ql, (jac, jac_t, jac_tt), U, V;
    constant_mass
  )
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  tfeop = TransientSemilinearFEOperator(
    mass, res_ql, U, V;
    constant_mass, order
  )
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  # TransientIMEXFEOperator
  constant_forms = (constant_mass,)

  im_tfeop = TransientSemilinearFEOperator(
    mass, im_res, (im_jac, im_jac_t, im_jac_tt), U, V;
    constant_mass
  )
  ex_tfeop = TransientFEOperator(ex_res, (ex_jac, ex_jac_t), U, V)
  tfeop = TransientIMEXFEOperator(im_tfeop, ex_tfeop)
  test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

  for constant_damping in (true, false)
    for constant_stiffness in (true, false)
      # TransientLinearFEOperator
      constant_forms = (constant_stiffness, constant_damping, constant_mass)

      tfeop = TransientLinearFEOperator(
        (stiffness, damping, mass), res_l, (jac, jac_t, jac_tt), U, V;
        constant_forms
      )
      test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)

      tfeop = TransientLinearFEOperator(
        (stiffness, damping, mass), res_l, U, V;
        constant_forms
      )
      test_transient_operator(tfeop, uhs0_fe, uhs0_cf, uhs0_dof, constant_forms)
    end
  end
end

end # module TransientFEOperatorsSolutionsTests
