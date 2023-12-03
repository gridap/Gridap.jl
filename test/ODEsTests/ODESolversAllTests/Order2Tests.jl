module Order2Tests

using Random
Random.seed!(123)
using Test

using Gridap
using Gridap.ODEs

include("../ODEOperatorsMocks.jl")
include("../ODESolversMocks.jl")

t0 = 0.0
dt = 1.0e-3
tF = t0 + 10 * dt

num_eqs = 2

M = randn(num_eqs, num_eqs)
while iszero(det(M))
  M = randn(num_eqs, num_eqs)
end
λ = randn(num_eqs)
μ = randn(num_eqs)
C = M * diagm(-(λ .+ μ))
K = M * diagm(λ .* μ)
α = randn(num_eqs)
f(t) = -M * exp.(α .* t)

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

odeop_nonlinear = ODEOperatorMock2{NonlinearODE}(M, C, K, f)

odeop_masslinear = ODEOperatorMock2{MassLinearODE}(M, C, K, f)

odeop_linear = ODEOperatorMock2{LinearODE}(M, C, K, f)
ODEs.is_jacobian_constant(odeop::typeof(odeop_linear), k::Integer) = true

odeops = [
  odeop_nonlinear,
  odeop_masslinear,
  odeop_linear
]

function test_solver(odeslvr, odeop, us0, tol)
  odesltn = solve(odeslvr, odeop, us0, t0, tF)

  for (uh_n, t_n) in odesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(abs2, eh_n))
    @test e_n < tol
  end
end

tol = 1.0e-4
disslvr_l = LUSolver()
disslvr_nl = NewtonRaphsonSolver(disslvr_l, 1.0e-8, 100)

# Solvers without memory
odeslvrs = [
]

us0 = (u0, v0,)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

# Solvers with memory
odeslvrs = [
  GeneralizedAlpha2(disslvr_nl, dt, 0.0),
  GeneralizedAlpha2(disslvr_nl, dt, 0.5),
  GeneralizedAlpha2(disslvr_nl, dt, 1.0),
  Newmark(disslvr_nl, dt, 0.5, 0.0),
  Newmark(disslvr_nl, dt, 0.5, 0.25),
]

a0 = -M \ (C * v0 + K * u0 + f(t0))
us0 = (u0, v0, a0)
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, us0, tol)
  end
end

end # module Order2Tests
