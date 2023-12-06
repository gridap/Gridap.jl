module Order1ODETests

using Test

using Gridap
using Gridap.ODEs

include("../ODEOperatorsMocks.jl")
include("../ODESolversMocks.jl")

t0 = 0.0
dt = 1.0e-3
tF = t0 + 10 * dt

num_eqs = 5

M = randn(num_eqs, num_eqs)
λ = randn(num_eqs)
K = M * diagm(-λ)

mass(t) = M
stiffness(t) = K
forms = (stiffness, mass)
form_zero(t) = zeros(num_eqs, num_eqs)

α = randn(num_eqs)
forcing(t) = -M * exp.(α .* t)
forcing_zero(t) = zero(t)

u0 = randn(num_eqs)

function u(t)
  s = zeros(typeof(t), num_eqs)
  for i in 1:num_eqs
    # Homogeneous solution
    s[i] += exp(λ[i] * t) * u0[i]
    # Particular solution
    s[i] += (exp(λ[i] * t) - exp(α[i] * t)) / (λ[i] - α[i])
  end
  s
end

odeop_nl = ODEOperatorMock{NonlinearODE}(forms, forcing)
odeop_ql = ODEOperatorMock{QuasilinearODE}(forms, forcing)
odeop_sl = ODEOperatorMock{SemilinearODE}(forms, forcing)
odeop_l = ODEOperatorMock{LinearODE}(forms, forcing)

# Testing some random combinations of `ODEOperatorType`s
odeop_imex1 = GenericIMEXODEOperator(
  ODEOperatorMock{LinearODE}((form_zero, mass), forcing_zero),
  ODEOperatorMock{QuasilinearODE}((stiffness,), forcing)
)

odeop_imex2 = GenericIMEXODEOperator(
  ODEOperatorMock{SemilinearODE}((form_zero, mass), forcing_zero),
  ODEOperatorMock{NonlinearODE}((stiffness,), forcing)
)

odeops = (
  odeop_nl,
  odeop_ql,
  odeop_sl,
  odeop_l,
  odeop_imex1,
  odeop_imex2,
)

function test_solver(odeslvr, odeop, us0, tol)
  odesltn = solve(odeslvr, odeop, us0, t0, tF)

  for (uh_n, t_n) in odesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(abs2, eh_n))
    @test e_n < tol
  end
end

tol = 1.0e-4
atol = 1.0e-12
rtol = 1.0e-8
maxiter = 100
disslvr_l = LUSolver()
disslvr_nl = DiscreteODESolverMock(rtol, atol, maxiter)

# Solvers without memory
odeslvrs = (
  ForwardEuler(disslvr_nl, dt),
  ThetaMethod(disslvr_nl, dt, 0.2),
  MidPoint(disslvr_nl, dt),
  ThetaMethod(disslvr_nl, dt, 0.8),
  BackwardEuler(disslvr_nl, dt),
  RungeKutta(disslvr_nl, disslvr_l, dt, :FE_1_0_1),
  RungeKutta(disslvr_nl, disslvr_l, dt, :SSPRK_3_0_3),
  RungeKutta(disslvr_nl, disslvr_l, dt, :BE_1_0_1),
  RungeKutta(disslvr_nl, disslvr_l, dt, :CN_2_0_2),
  RungeKutta(disslvr_nl, disslvr_l, dt, :SDIRK_2_0_2),
  RungeKutta(disslvr_nl, disslvr_l, dt, :SDIRK_2_0_3),
  RungeKutta(disslvr_nl, disslvr_l, dt, :ESDIRK_3_1_2),
  RungeKutta(disslvr_nl, disslvr_l, dt, :TRBDF2_3_2_3),
  RungeKutta(disslvr_nl, disslvr_l, dt, :TRX2_3_2_3),
)

us0 = (u0,)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

# Solvers with memory
odeslvrs = [
  GeneralizedAlpha1(disslvr_nl, dt, 0.0),
  GeneralizedAlpha1(disslvr_nl, dt, 0.5),
  GeneralizedAlpha1(disslvr_nl, dt, 1.0),
]

v0 = -M \ (K * u0 + forcing(t0))
us0 = (u0, v0,)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

# Solvers for `IMEXODEOperator`s
odeops = (
  odeop_imex1,
  odeop_imex2,
)

odeslvrs = (
  RungeKutta(disslvr_nl, disslvr_l, dt, :IMEX_FE_BE_2_0_1),
  RungeKutta(disslvr_nl, disslvr_l, dt, :IMEX_Midpoint_2_0_2),
)

us0 = (u0,)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

end # module Order1ODETests
