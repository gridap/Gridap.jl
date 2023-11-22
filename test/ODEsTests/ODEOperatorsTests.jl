module ODEOperatorsTests

using Test

using SparseArrays: spzeros

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")

a, b, c = 1.0, 2.0, 3.0
op = ODEOperatorMock{ConstantMassODE}(a, b, c, 1)

t = 0.0
u = ones(2)
u̇ = 2 * ones(2)

cache = allocate_cache(op)
update_cache!(cache, op, t)

r = allocate_residual(op, t, (u, u̇), cache)
@test size(r) == (2,)

J = allocate_jacobian(op, t, (u, u̇), cache)
@test size(J) == (2, 2)

residual!(r, op, t, (u, u̇), cache)
_r = zeros(2)
_r[1] = u̇[1] - a * u[1]
_r[2] = u̇[2] - b * u[1] - c * u[2]
@test all(r .== _r)

fill!(J, 0)
jacobian!(J, op, t, (u, u̇), 0, 1, cache)
_J = zeros(2, 2)
_J[1, 1] = -a
_J[2, 1] = -b
_J[2, 2] = -c
@test all(J .== _J)

jacobian!(J, op, t, (u, u̇), 1, 1, cache)
_J[1, 1] += 1
_J[2, 2] += 1
@test all(J .== _J)

@test test_ode_operator(op, t, (u, u̇))

end # module ODEOperatorsTests
