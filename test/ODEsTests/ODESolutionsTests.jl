module ODESolutionsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

t0 = randn()
tF = t0 + rand()
dt = (tF - t0) / 10
u0 = randn(2)

M = randn(2, 2)
K = randn(2, 2)
while iszero(det(M + dt * K))
  M = randn(2, 2)
  K = randn(2, 2)
end
f(t) = [cospi(t), sinpi(t)]
ode_op = ODEOperatorMock{MassLinearODE}(M, K, f)

nls = NLSolverMock()
solver = ODESolverMock(nls, dt)

sol = solve(solver, ode_op, u0, t0, tF)

_J = M + dt * K
tprev = t0
uprev = copy(u0)
for (u_n, t_n) in sol
  global tprev, uprev
  @test t_n ≈ tprev + dt
  @test u_n ≈ uprev + dt * (-_J \ (K * uprev + f(t_n)))
  tprev = t_n
  copy!(uprev, u_n)
end

@test test_ode_solution(sol)

end # module ODESolutionsTests
