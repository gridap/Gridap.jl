module NLSolversTests

using Gridap.Algebra
using Gridap.Algebra: NonLinearOperatorMock

op = NonLinearOperatorMock()
nls = NLSolver(show_trace=false,method=:newton)

x0 = zero_initial_guess(op)
x = [1.0, 3.0]
test_non_linear_solver(nls,op,x0,x)

x0 = [2.1,2.9]
x = [2.0, 3.0]
test_non_linear_solver(nls,op,x0,x)

x0 = zero_initial_guess(op)
cache = solve!(x0,nls,op)

_ = cache.result

end # module
