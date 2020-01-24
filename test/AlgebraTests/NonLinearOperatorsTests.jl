module NonLinearOperatorsTests

using Test
using Gridap.Algebra
using Gridap.Algebra: NonLinearOperatorMock

op = NonLinearOperatorMock()

x = [1.0, 3.0]
b = [0.0,0.0]
A = [-1.0 0.0; 0.0 1.0]

test_non_linear_operator(op,x,b,jac=A)

@test residual(op,x) ≈ b
@test jacobian(op,x) ≈ A
b1, A1 = residual_and_jacobian(op,x)
@test b1 ≈ b
@test A1 ≈ A

end # module
