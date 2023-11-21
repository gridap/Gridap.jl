module ODESolutionsTests

using Test

using Gridap
using Gridap.ODEs
using Gridap.ODEs.ODETools
using Gridap.ODEs.ODETools: GenericODESolution

include("ODEOperatorMocks.jl")
include("ODESolverMocks.jl")

a, b, c = 1.0, 0.0, 1.0
op = ODEOperatorMock{ConstantMassODE}(a, b, c, 1)

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

steps = solve(solver, op, u0, t0, tF)

uF = copy(u0)
@. uF = 1.0

for (n, (u_n, t_n)) in enumerate(steps)
  @test t_n ≈ t0 + n * dt
  if a == c
    D = α^n
    @test u_n[1] ≈ u0[1] / D
    @test u_n[2] ≈ (u0[2] - n * β * u0[1]) / D
  else
    D = (α * γ)^n
    @test u_n[1] ≈ (γ^n * u0[1]) / D
    @test u_n[2] ≈ (α^n * u0[2] - β * (α^n - γ^n) / (α - γ)) / D
  end
end

@test test_ode_solution(steps)

end # module ODESolutionsTests
