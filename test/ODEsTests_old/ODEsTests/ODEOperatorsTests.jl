module ODEOperatorsTests

using Gridap.ODEs.ODETools
using Test

import Gridap.ODEs.ODETools: test_ode_operator

include("ODEOperatorMocks.jl")

op = ODEOperatorMock{Float64,Constant}(1.0,2.0,3.0,1)

u = ones(2)
u_t = ones(2)*2.0

@assert(length(u) == 2)
@assert(length(u_t) == 2)

cache = allocate_cache(op)
update_cache!(cache,op,0.0)

t = 0.0
r = allocate_residual(op,t,u,cache)
@test r == zeros(2)

J = allocate_jacobian(op,t,u,cache)
@test J == zeros(2,2)

residual!(r,op,t,(u,u_t),cache)
_r = zeros(2)
_r[1] = u_t[1] - op.a * u[1]
_r[2] = u_t[2] - op.b * u[1] - op.c * u[2]
@test all(r .== _r)

J .= 0
jacobian!(J,op,t,(u,u_t),1,1.0,cache)
_J = zeros(2,2)
_J[1,1] = -op.a
_J[2,1] = -op.b
_J[2,2] = -op.c
@test all(J .== _J)

jacobian!(J,op,t,(u,u_t),2,1.0,cache)
_J[1,1] += 1.0
_J[2,2] += 1.0
@test all(J .== _J)
_J
J

@test test_ode_operator(op,t,u,u_t)

end #module
