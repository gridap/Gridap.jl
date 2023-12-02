module ODESolutionsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

t0 = randn()
tF = t0 + rand()
dt = (tF - t0) / 10

num_eqs = 2

M = randn(num_eqs, num_eqs)
C = randn(num_eqs, num_eqs)
K = randn(num_eqs, num_eqs)
while iszero(det(M + dt * K)) || iszero(det(M + dt * C + 3 * dt^2 / 2 * K))
  M = randn(num_eqs, num_eqs)
  C = randn(num_eqs, num_eqs)
  K = randn(num_eqs, num_eqs)
end
α = randn(num_eqs)
f(t) = -exp.(-α .* t)

u0 = randn(num_eqs)
v0 = randn(num_eqs)

_odeop1 = ODEOperatorMock1
_odeslvr1 = ODESolverMock1
forms1 = (M, K)
us01 = (u0,)
function _u1(t, (u,))
  x = -(M + dt * K) \ (K * u + f(t + dt))
  u + dt * x
end

_odeop2 = ODEOperatorMock2
_odeslvr2 = ODESolverMock2
forms2 = (M, C, K)
us02 = (u0, v0,)
function _u2(t, (u, v))
  x = -(M + dt * C + 3 * dt^2 / 2 * K) \ (C * v + K * (u + dt * v) + f(t + dt))
  u + dt * v + (3 * dt^2 / 2) * x
end

disslvr = NLSolverMock()

for (_odeop, _odeslvr, forms, us0, _u) in (
  (_odeop1, _odeslvr1, forms1, us01, _u1),
  (_odeop2, _odeslvr2, forms2, us02, _u2),
)
  for T in (NonlinearODE, MassLinearODE, LinearODE,)
    odeop = _odeop{T}(forms..., f)
    odeslvr = _odeslvr(disslvr, dt)

    odesltn = solve(odeslvr, odeop, us0, t0, tF)

    tprev = t0
    uprev = copy.(us0)

    it = iterate(odesltn)
    while !isnothing(it)
      data, state = it
      u_n, t_n = data
      us_n, = state

      @test t_n ≈ tprev + dt
      @test u_n ≈ _u(tprev, uprev)

      tprev = t_n
      for (uprevi, us_ni) in zip(uprev, us_n)
        copy!(uprevi, us_ni)
      end

      it = iterate(odesltn, state)
    end

    @test test_ode_solution(odesltn)
  end
end

end # module ODESolutionsTests
