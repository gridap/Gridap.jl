module Order1ODETests

using Test
using SparseArrays

using Gridap
using Gridap.ODEs

include("../ODEOperatorsMocks.jl")
include("../ODESolversMocks.jl")

# M u̇ + M * Diag(-λ) u = M * f(t),
# where f(t) = exp(Diag(α) * t)
t0 = 0.0
dt = 1.0e-3
tF = t0 + 10 * dt

num_eqs = 5

M = sprandn(num_eqs, num_eqs, 1.0)
λ = randn(num_eqs)
K = M * spdiagm(-λ)

mass(t) = M
stiffness(t) = K
forms = (stiffness, mass)
mat0 = sprand(num_eqs, num_eqs, 1.0)
nonzeros(mat0) .= 0
form_zero(t) = mat0

α = randn(num_eqs)
forcing(t) = M * exp.(α .* t)
forcing_zero(t) = zeros(typeof(t), num_eqs)

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
  odesltn = solve(odeslvr, odeop, t0, tF, us0)

  for (t_n, uh_n) in odesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(abs2, eh_n))
    @test e_n < tol
  end
end

tol = 1.0e-4
atol = 1.0e-12
rtol = 1.0e-8
maxiter = 100
sysslvr_l = LUSolver()
sysslvr_nl = NonlinearSolverMock(rtol, atol, maxiter)

odeslvrs = (
  ForwardEuler(sysslvr_nl, dt),
  ThetaMethod(sysslvr_nl, dt, 0.2),
  MidPoint(sysslvr_nl, dt),
  ThetaMethod(sysslvr_nl, dt, 0.8),
  BackwardEuler(sysslvr_nl, dt),
  GeneralizedAlpha1(sysslvr_nl, dt, 0.0),
  GeneralizedAlpha1(sysslvr_nl, dt, 0.5),
  GeneralizedAlpha1(sysslvr_nl, dt, 1.0),
)
for tableau in available_tableaus
  global odeslvrs
  odeslvr = RungeKutta(sysslvr_nl, sysslvr_l, dt, tableau)
  odeslvrs = (odeslvrs..., odeslvr)
end

us0 = (u0,)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

# Tests with initial velocity
odeslvrs = (
  GeneralizedAlpha1(sysslvr_nl, dt, 0.0),
  GeneralizedAlpha1(sysslvr_nl, dt, 0.5),
  GeneralizedAlpha1(sysslvr_nl, dt, 1.0),
)

v0 = exp.(α .* t0) - spdiagm(-λ) * u0
us0 = (u0, v0)
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

odeslvrs = ()
for tableau in available_imex_tableaus
  global odeslvrs
  odeslvr = RungeKutta(sysslvr_nl, sysslvr_l, dt, tableau)
  odeslvrs = (odeslvrs..., odeslvr)
end

us0 = (u0,)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

end # module Order1ODETests
