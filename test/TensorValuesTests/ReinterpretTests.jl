module ReinterpretTests

using Test
using Gridap.TensorValues

v = VectorValue(1.3,2.1)

V = fill(v,5)

A = reinterpret(V)

R = [1.3 1.3 1.3 1.3 1.3; 2.1 2.1 2.1 2.1 2.1]

@test A == R

end # module ReinterpretTests
