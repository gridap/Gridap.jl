module Order2ODETests

using Test
using SparseArrays

using Gridap
using Gridap.ODEs

include("../ODEOperatorsMocks.jl")
include("../ODESolversMocks.jl")

# M ü + M * Diag(-(λ .+ μ)) u̇ + M * Diag(λ .* μ) u = M * f(t),
# where f(t) = exp(Diag(α) * t)
t0 = 0.0
dt = 1.0e-3
tF = t0 + 10 * dt

num_eqs = 5

M = sprandn(num_eqs, num_eqs, 1.0)
λ = randn(num_eqs)
μ = randn(num_eqs)
C = M * spdiagm(-(λ .+ μ))
K = M * spdiagm(λ .* μ)

mass(t) = M
damping(t) = C
stiffness(t) = K
forms = (stiffness, damping, mass)
mat0 = sprand(num_eqs, num_eqs, 1.0)
nonzeros(mat0) .= 0
form_zero(t) = mat0

α = randn(num_eqs)
forcing(t) = M * exp.(α .* t)
forcing_zero(t) = zeros(typeof(t), num_eqs)

u0 = randn(num_eqs)
v0 = randn(num_eqs)

function u(t)
  s = zeros(typeof(t), num_eqs)
  for i in 1:num_eqs
    # Homogeneous solution
    s[i] += (μ[i] * u0[i] - v0[i]) / (μ[i] - λ[i]) * exp(λ[i] * t)
    s[i] += (λ[i] * u0[i] - v0[i]) / (λ[i] - μ[i]) * exp(μ[i] * t)
    # Particular solution
    s[i] += (exp(λ[i] * t) - exp(α[i] * t)) / (λ[i] - μ[i]) / (λ[i] - α[i])
    s[i] += (exp(μ[i] * t) - exp(α[i] * t)) / (μ[i] - λ[i]) / (μ[i] - α[i])
  end
  s
end

odeop_nl = ODEOperatorMock{NonlinearODE}(forms, forcing)
odeop_ql = ODEOperatorMock{QuasilinearODE}(forms, forcing)
odeop_sl = ODEOperatorMock{SemilinearODE}(forms, forcing)
odeop_l = ODEOperatorMock{LinearODE}(forms, forcing)

# Testing some random combinations of `IMEXODEOperator`s
odeop_imex1 = GenericIMEXODEOperator(
  ODEOperatorMock{LinearODE}((stiffness, form_zero, mass), forcing_zero),
  ODEOperatorMock{QuasilinearODE}((form_zero, damping,), forcing)
)

odeop_imex2 = GenericIMEXODEOperator(
  ODEOperatorMock{SemilinearODE}((form_zero, damping, mass), forcing_zero),
  ODEOperatorMock{NonlinearODE}((stiffness, form_zero,), forcing)
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
  GeneralizedAlpha2(sysslvr_nl, dt, 0.0),
  GeneralizedAlpha2(sysslvr_nl, dt, 0.5),
  GeneralizedAlpha2(sysslvr_nl, dt, 1.0),
  Newmark(sysslvr_nl, dt, 0.5, 0.0),
  Newmark(sysslvr_nl, dt, 0.5, 0.25),
)

us0 = (u0, v0)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

# Tests with initial acceleration
a0 = exp.(α .* t0) - spdiagm(-(λ .+ μ)) * v0 - spdiagm(λ .* μ) * u0
us0 = (u0, v0, a0)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

end # module Order2ODETests
