module BubbleMonomialBasesTests

using Test
using Gridap.Polynomials
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: Broadcasting

# 2D Bubble Monomial Basis
x1 = Point(0.0, 0.0)
x2 = Point(1.0, 0.0)
x3 = Point(0.0, 1.0)
x4 = Point(0.1, 0.1)
x0 = sum([x1, x2, x3]) / 3.0
g((x, y)) = x * y * (1.0 - x - y) * 27.0
∇g((x, y)) = VectorValue((y - 2.0*x*y-y^2)*27.0, (x - 2.0*x*y-x^2)*27.0)

T = Float64
G = gradient_type(T, x0)
f = BubbleMonomialBasis(T, Val(2))
∇f = Broadcasting(∇)(f)
@test evaluate(f, [x0, x1, x2, x3]) ≈ map(g, [x0, x1, x2, x3])
@test evaluate(∇f, [x0, x4]) ≈ map(∇g, [x0; x4;;])

T = VectorValue{2, Float64}
G = gradient_type(T, x0)
f = BubbleMonomialBasis(T, Val(2))
∇f = Broadcasting(∇)(f)
@test evaluate(f, [x0, x1, x2, x3]) ≈ T[
	(1.0, 0.0) (0.0, 1.0);
	(0.0, 0.0) (0.0, 0.0);
	(0.0, 0.0) (0.0, 0.0);
	(0.0, 0.0) (0.0, 0.0);
]
@test evaluate(∇f, [x0, x4]) ≈ G[
	(0.0, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 0.0);
	(1.89, 1.89, 0.0, 0.0) (0.0, 0.0, 1.89, 1.89);
]


# 3D Bubble Monomial Basis
x1 = Point(0.0, 0.0, 0.0)
x2 = Point(1.0, 0.0, 0.0)
x3 = Point(0.0, 1.0, 0.0)
x4 = Point(0.0, 0.0, 1.0)
x0 = sum([x1, x2, x3, x4]) / 4.0

T = Float64
f = BubbleMonomialBasis(T, Val(3))
@test evaluate(f, [x0, x1, x2, x3, x4]) ≈ [1.0, 0.0, 0.0, 0.0, 0.0]

T = VectorValue{3, Float64}
f = BubbleMonomialBasis(T, Val(3))
@test evaluate(f, [x0, x1, x2, x3, x4]) ≈ T[
	(1.0, 0.0, 0.0) (0.0, 1.0, 0.0) (0.0, 0.0, 1.0);
	(0.0, 0.0, 0.0) (0.0, 0.0, 0.0) (0.0, 0.0, 0.0);
	(0.0, 0.0, 0.0) (0.0, 0.0, 0.0) (0.0, 0.0, 0.0);
	(0.0, 0.0, 0.0) (0.0, 0.0, 0.0) (0.0, 0.0, 0.0);
	(0.0, 0.0, 0.0) (0.0, 0.0, 0.0) (0.0, 0.0, 0.0);
]

end
