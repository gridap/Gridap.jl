module ODEOperatorsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")

M = randn(2, 2)
K = randn(2, 2)
f(t) = [cospi(t), sinpi(t)]
op = ODEOperatorMock{MassLinearODE}(M, K, f)

t = randn()
u = randn(2)
u̇ = randn(2)

cache = allocate_cache(op, t, (u, u̇))
update_cache!(cache, op, t)

r = allocate_residual(op, t, (u, u̇), cache)
@test size(r) == (2,)

J = allocate_jacobian(op, t, (u, u̇), cache)
@test size(J) == (2, 2)

_r = M * u̇ + K * u + f(t)
residual!(r, op, t, (u, u̇), cache)
@test r ≈ _r

_J = K
fill!(J, 0)
jacobian!(J, op, t, (u, u̇), 0, 1, cache)
@test J ≈ _J

_J += M
jacobian!(J, op, t, (u, u̇), 1, 1, cache)
@test J ≈ _J

@test test_ode_operator(op, t, (u, u̇))

end # module ODEOperatorsTests
