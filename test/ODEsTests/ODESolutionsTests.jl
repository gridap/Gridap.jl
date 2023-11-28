module ODESolutionsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

a, b, c = 1.0, 0.0, 1.0
op = ODEOperatorMock{MassLinearODE}(a, b, c, 1)

t0 = 0.0
tF = 1.0
dt = 0.1
u0 = 2 * ones(2)

α = 1 - dt * a
β = -dt * b
γ = 1 - dt * c
@assert !iszero(α)
@assert !iszero(γ)

nls = NLSolverMock()
solver = ODESolverMock(nls, dt)

sol = solve(solver, op, u0, t0, tF)

uF = copy(u0)
fill!(uF, 1)

αⁿ, γⁿ = one(α), one(γ)
for (n, (u_n, t_n)) in enumerate(sol)
  @test t_n ≈ t0 + n * dt
  if a == c
    D = α^n
    @test u_n[1] ≈ u0[1] / D
    @test u_n[2] ≈ (u0[2] - n * β * u0[1]) / D
  else
    D = αⁿ * γⁿ
    αⁿ = α * αⁿ
    γⁿ = γ * γⁿ
    @test u_n[1] ≈ (γⁿ * u0[1]) / D
    @test u_n[2] ≈ (αⁿ * u0[2] - β * (αⁿ - γⁿ) / (α - γ)) / D
  end
end

@test test_ode_solution(sol)

end # module ODESolutionsTests
