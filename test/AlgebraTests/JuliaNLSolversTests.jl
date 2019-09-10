module JuliaNLSolversTests

using Test
using Gridap

using ..NonLinearOperatorMocks

op = NonLinearOperatorMock()

nls = JuliaNLSolver(show_trace=false,method=:newton)

x = solve(nls,op)
@test x ≈ [1.0, 3.0]

cache = solve!(x,nls,op)

solve!(x,nls,op,cache)

cache.result

end # module
