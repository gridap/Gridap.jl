module NonLinearSolversTests

using Gridap.Algebra
using Gridap.Algebra: NonLinearOperatorMock

op = NonLinearOperatorMock()

ls = LUSolver()
tol = 1.e-10
maxiters = 20
nls = NewtonRaphsonSolver(ls,tol,maxiters)

x0 = zero_initial_guess(op)
x = [1.0, 3.0]
test_non_linear_solver(nls,op,x0,x)

x0 = [2.1,2.9]
x = [2.0, 3.0]
test_non_linear_solver(nls,op,x0,x)

end # module
