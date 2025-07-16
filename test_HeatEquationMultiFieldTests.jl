# module HeatEquationMultifieldTests

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
V = FESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)

Y = MultiFieldFESpace([V, V])
X = TransientMultiFieldFESpace([U, U])

# Integration
Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

# FE operator
ft(t) = x -> ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)
_mass(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
_mass(t, u, ∂ₜu, v) = _mass(t, ∂ₜu, v)
_stiffness(t, u, v) = ∫(∇(u) ⊙ ∇(v)) * dΩ
_forcing(t, v) = ∫(f(t) ⋅ v) * dΩ

_res(t, u, v) = _mass(t, ∂t(u), v) + _stiffness(t, u, v) - _forcing(t, v)
_res_ql(t, u, v) = _stiffness(t, u, v) - _forcing(t, v)
_res_l(t, v) = _forcing(t, v)

mass(t, (∂ₜu1, ∂ₜu2), (v1, v2)) = _mass(t, ∂ₜu1, v1) + _mass(t, ∂ₜu2, v2)
mass(t, (u1, u2), (∂ₜu1, ∂ₜu2), (v1, v2)) = _mass(t, u1, ∂ₜu1, v1) + _mass(t, u2, ∂ₜu2, v2)
stiffness(t, (u1, u2), (v1, v2)) = _stiffness(t, u1, v1) + _stiffness(t, u2, v2)

res(t, (u1, u2), (v1, v2)) = _res(t, u1, v1) + _res(t, u2, v2)
jac(t, x, (du1, du2), (v1, v2)) = _stiffness(t, du1, v1) + _stiffness(t, du2, v2)
jac_t(t, x, (dut1, dut2), (v1, v2)) = _mass(t, dut1, v1) + _mass(t, dut2, v2)

res_ql(t, (u1, u2), (v1, v2)) = _res_ql(t, u1, v1) + _res_ql(t, u2, v2)
res_l(t, (v1, v2)) = _res_l(t, v1) + _res_l(t, v2)

args_ad = (X, Y)
tfeop_nl_ad = TransientFEOperator(res, args_ad...)

using Gridap.ODEs: get_order, _make_uh_from_us, _to_transient_single_fields
using Gridap.CellData

odeop = ODEOpFromTFEOp(tfeop_nl_ad)
order = get_order(odeop)
Ut = get_trial(odeop.tfeop)
U = allocate_space(Ut)
us = (get_free_dof_values(zero(U(0.0))),get_free_dof_values(zero(U(0.0))))
Uts = (Ut,)
Us = (U,)
for k in 1:order
  Uts = (Uts..., ∂t(Uts[k]))
  Us = (Us..., allocate_space(Uts[k+1]))
end

# Variables for assembly
uh = _make_uh_from_us(odeop, us, Us)
V = get_test(odeop.tfeop)
v = get_fe_basis(V)
Ut = evaluate(get_trial(odeop.tfeop), nothing)
du = get_trial_fe_basis(Ut)
jacs = get_jacs(odeop.tfeop)

# function ODEs.TransientCellField(cell_fields::Tuple, derivatives::Tuple)
#   _to_transient_single_fields(cell_fields,derivatives)
# end

# _tst = (uh.derivatives[1],uh.derivatives[1])

jacs[1](0.0, uh, du, v)

# TransientCellField((uh.cellfield[1],uh.cellfield[2]),uh.derivatives, U)

# _to_transient_single_fields((uh.cellfield[1],uh.cellfield[2]),uh.derivatives)

nothing







tfeops = (
  tfeop_nl_man,
  tfeop_ql_man,
  tfeop_sl_man,
  tfeop_l_man,
  tfeop_nl_ad,
  tfeop_ql_ad,
  tfeop_sl_ad,
  tfeop_l_ad,
)

# Initial conditions
t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(t0)
uh0 = interpolate_everywhere(u(t0), U0)
X0 = X(t0)
xh0 = interpolate_everywhere([uh0, uh0], X0)
xhs0 = (xh0,)

_tst = ODEOpFromTFEOp(tfeop_nl_ad)
opcache=allocate_odeopcache(ODEOpFromTFEOp(tfeop_nl_man),0.0,(get_free_dof_values(xh0),get_free_dof_values(xh0)))
allocate_jacobian(ODEOpFromTFEOp(tfeop_nl_man),0.0,(get_free_dof_values(xh0),get_free_dof_values(xh0)),opcache)

allocate_odeopcache(_tst, 0.0, (get_free_dof_values(xh0),get_free_dof_values(xh0)))

allocate_odecache(odeslvrs[1],_tst,0.0,(get_free_dof_values(xh0),))

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
    fesltn = solve(odeslvr, tfeop, t0, tF, xhs0)

    for (t_n, xhs_n) in fesltn
      eh_n = u(t_n) - xhs_n[1]
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol

      eh_n = u(t_n) - xhs_n[2]
      e_n = sqrt(sum(∫(eh_n ⋅ eh_n) * dΩ))
      @test e_n < tol
    end
  end
end

# end # module HeatEquationMultifieldTests
