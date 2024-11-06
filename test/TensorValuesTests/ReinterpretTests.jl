module ReinterpretTests

using Test
using Gridap.TensorValues

v = VectorValue(1.3,2.1)
V = fill(v,5)
A = reinterpret(V)
R = [1.3 1.3 1.3 1.3 1.3; 2.1 2.1 2.1 2.1 2.1]
@test A == R

v = TensorValue(2.0,3.5,1.0,4.0)
V = fill(v,5)
A = reinterpret(V)
R = [2.0 2.0 2.0 2.0 2.0; 3.5 3.5 3.5 3.5 3.5; 1.0 1.0 1.0 1.0 1.0; 4.0 4.0 4.0 4.0 4.0]
@test A == R

v = SymFourthOrderTensorValue(1,2,3, 4,5,6, 7,8,9)
V = fill(v,5)
A = reinterpret(V)
R = [1 1 1 1 1; 2 2 2 2 2; 3 3 3 3 3; 4 4 4 4 4; 5 5 5 5 5; 6 6 6 6 6; 7 7 7 7 7; 8 8 8 8 8; 9 9 9 9 9]
@test A == R

v = SymTensorValue(1,2,3)
V = fill(v,5)
A = reinterpret(V)
R = [1 1 1 1 1; 2 2 2 2 2; 3 3 3 3 3]
@test A == R

v = SymTracelessTensorValue(1,2)
V = fill(v,5)
A = reinterpret(V)
R = [1 1 1 1 1; 2 2 2 2 2; -1 -1 -1 -1 -1]
@test A == R

end # module ReinterpretTests
