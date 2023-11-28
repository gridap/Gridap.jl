module TimeDerivativesTests

using Test

using ForwardDiff

using Gridap
using Gridap.ODEs

# First time derivative, scalar-valued
f1(x, t) = 5 * x[1] * x[2] + x[2]^2 * t^3
∂tf1(x, t) = 3 * x[2]^2 * t^2

f2(x, t) = t^2
∂tf2(x, t) = 2 * t

f3(x, t) = x[1]^2
∂tf3(x, t) = zero(x[1])

f4(x, t) = x[1]^t^2
∂tf4(x, t) = 2 * t * log(x[1]) * f4(x, t)

for (f, ∂tf) in ((f1, ∂tf1), (f2, ∂tf2), (f3, ∂tf3), (f4, ∂tf4),)
  dtf = (x, t) -> ForwardDiff.derivative(t -> f(x, t), t)

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂t(f)(xv, tv) ≈ ∂tf(xv, tv)
  @test ∂t(f)(xv, tv) ≈ dtf(xv, tv)
  @test ∂t(f)(xv, tv) ≈ ∂t(f)(xv)(tv) ≈ ∂t(f)(tv)(xv)
end

# First time derivative, vector-valued
f1(x, t) = VectorValue(5 * x[1] * x[2], x[2]^2 * t^3)
∂tf1(x, t) = VectorValue(zero(x[1]), x[2]^2 * 3 * t^2)

f2(x, t) = VectorValue(x[1]^2, zero(x[2]))
∂tf2(x, t) = VectorValue(zero(x[1]), zero(x[2]))

f3(x, t) = VectorValue(x[1]^2, t)
∂tf3(x, t) = VectorValue(zero(x[1]), one(t))

for (f, ∂tf) in ((f1, ∂tf1), (f2, ∂tf2), (f3, ∂tf3),)
  dtf = (x, t) -> VectorValue(ForwardDiff.derivative(t -> get_array(f(x, t)), t))

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂t(f)(xv, tv) ≈ ∂tf(xv, tv)
  @test ∂t(f)(xv, tv) ≈ dtf(xv, tv)
  @test ∂t(f)(xv, tv) ≈ ∂t(f)(xv)(tv) ≈ ∂t(f)(tv)(xv)
end

# First time derivative, tensor-valued
f1(x, t) = TensorValue(x[1] * t, x[2] * t, x[1] * x[2], x[1] * t^2)
∂tf1(x, t) = TensorValue(x[1], x[2], zero(x[1]), 2 * x[1] * t)

for (f, ∂tf) in ((f1, ∂tf1),)
  dtf = (x, t) -> TensorValue(ForwardDiff.derivative(t -> get_array(f(x, t)), t))

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂t(f)(xv, tv) ≈ ∂tf(xv, tv)
  @test ∂t(f)(xv, tv) ≈ dtf(xv, tv)
  @test ∂t(f)(xv, tv) ≈ ∂t(f)(xv)(tv) ≈ ∂t(f)(tv)(xv)
end

# Spatial derivatives
# f(x, t) = VectorValue(x[1]^2, t)
# ∇f(x, t) = ∇(y -> f(y, t))(x)

# tv = rand(Float64)
# xv = Point(rand(Float64, 2)...)
# ∇f(xv, tv)

# TODO
# @santiagobadia : Is there any way to make this transparent to the user
# I guess not unless we create a type for these analytical (space-only or
# space-time via a trait) functions
# Probably a try-catch?

# Second time derivative, scalar-valued
f1(x, t) = t^2
∂tf1(x, t) = 2 * t
∂ttf1(x, t) = 2 * one(t)

f2(x, t) = x[1] * t^2
∂tf2(x, t) = 2 * x[1] * t
∂ttf2(x, t) = 2 * x[1]

for (f, ∂tf, ∂ttf) in ((f1, ∂tf1, ∂ttf1), (f2, ∂tf2, ∂ttf2),)
  dtf = (x, t) -> ForwardDiff.derivative(t -> f(x, t), t)
  dttf = (x, t) -> ForwardDiff.derivative(t -> dtf(x, t), t)

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂tt(f)(xv, tv) ≈ ∂ttf(xv, tv)
  @test ∂tt(f)(xv, tv) ≈ dttf(xv, tv)
  @test ∂tt(f)(xv, tv) ≈ ∂tt(f)(xv)(tv) ≈ ∂tt(f)(tv)(xv)
end

end # module TimeDerivativesTests
