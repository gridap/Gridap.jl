module NLSolversTests

using Test
using Gridap

using ..NonLinearOperatorMocks

op = NonLinearOperatorMock()

nls = NLSolver(show_trace=false,method=:newton)

x = solve(nls,op)
@test x ≈ [1.0, 3.0]

cache = solve!(x,nls,op)

@test string(cache) == "NLSolversCache object"

solve!(x,nls,op,cache)

cache.result

end # module
