module NonLinearSolversTests

using Test
using Gridap
using ..NonLinearOperatorMocks

op = NonLinearOperatorMock()

x = [2.0, 3.0]
x = [1.0, 3.0]
@test residual(op,x) ≈ [0.0,0.0]
@test jacobian(op,x) ≈ [-1.0 0.0; 0.0 1.0]
@test create_in_domain(op) == zeros(2)

ls = LUSolver()
tol = 1.e-10
maxiters = 20
nls = NewtonRaphsonSolver(ls,tol,maxiters)

x = solve(nls,op)
@test x ≈ [1.0, 3.0]

x = [2.1,2.9]
cache = solve!(x,nls,op)
@test x ≈ [2.0, 3.0]

x = [2.1,2.9]
solve!(x,nls,op,cache)
@test x ≈ [2.0, 3.0]

end
