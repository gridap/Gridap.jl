module RungeKuttaTests

using Test

using Gridap
using Gridap.ODEs

include("../ODEOperatorsMocks.jl")
include("../ODESolversMocks.jl")

t0 = 0.0
dt = 1.0e-3
tF = t0 + 10 * dt
u0 = randn(2)

M = randn(2, 2)
while iszero(det(M))
  M = randn(2, 2)
end
α, β = randn(), randn()
K = M * diagm([α, β])
f(t) = M * [cospi(t), sinpi(t)]

function u(t)
  s = zeros(typeof(t), 2)
  s[1] = exp(-α * t) * (u0[1] - (exp(α * t) * (α * cospi(t) + pi * sinpi(t)) - exp(α * t0) * (α * cospi(t0) + pi * sinpi(t0))) / (α^2 + π^2))
  s[2] = exp(-β * t) * (u0[2] - (exp(β * t) * (β * sinpi(t) - pi * cospi(t)) - exp(β * t0) * (β * sinpi(t0) - pi * cospi(t0))) / (β^2 + π^2))
  s
end

odeop_nonlinear = ODEOperatorMock1{NonlinearODE}(M, K, f)

odeop_masslinear = ODEOperatorMock1{MassLinearODE}(M, K, f)

odeop_linear = ODEOperatorMock1{LinearODE}(M, K, f)
ODEs.is_jacobian_constant(odeop::typeof(odeop_linear), k::Integer) = true

odeops = [
  odeop_nonlinear,
  odeop_masslinear,
  odeop_linear
]

function test_solver(odeslvr, odeop, tol)
  odesltn = solve(odeslvr, odeop, u0, t0, tF)

  for (uh_n, t_n) in odesltn
    eh_n = u(t_n) - uh_n
    e_n = sqrt(sum(abs2, eh_n))
    @test e_n < tol
  end
end

tol = 1.0e-4
disslvr_l = LUSolver()
disslvr_nl = NewtonRaphsonSolver(disslvr_l, 1.0e-8, 100)

odeslvrs = [
  RungeKutta(disslvr_nl, disslvr_l, dt, :FE_1_0_1)
  RungeKutta(disslvr_nl, disslvr_l, dt, :SSPRK_3_0_3)
  RungeKutta(disslvr_nl, disslvr_l, dt, :BE_1_0_1)
  RungeKutta(disslvr_nl, disslvr_l, dt, :CN_2_0_2)
  RungeKutta(disslvr_nl, disslvr_l, dt, :SDIRK_2_0_2)
  RungeKutta(disslvr_nl, disslvr_l, dt, :SDIRK_2_0_3)
  RungeKutta(disslvr_nl, disslvr_l, dt, :ESDIRK_3_1_2)
  RungeKutta(disslvr_nl, disslvr_l, dt, :TRBDF2_3_2_3)
  RungeKutta(disslvr_nl, disslvr_l, dt, :TRX2_3_2_3)
]

# Main loop
for odeslvr in odeslvrs
  for odeop in odeops
    test_solver(odeslvr, odeop, tol)
  end
end

end # module RungeKuttaTests
