module TimeDerivativesTests

using Test

using ForwardDiff

using Gridap
using Gridap.ODEs

# First time derivative, scalar-valued
f1(t) = x -> 5 * x[1] * x[2] + x[2]^2 * t^3
∂tf1(t) = x -> 3 * x[2]^2 * t^2

f2(t) = x -> t^2
∂tf2(t) = x -> 2 * t

f3(t) = x -> x[1]^2
∂tf3(t) = x -> zero(x[1])

f4(t) = x -> x[1]^t^2
∂tf4(t) = x -> 2 * t * log(x[1]) * f4(t)(x)

for (f, ∂tf) in ((f1, ∂tf1), (f2, ∂tf2), (f3, ∂tf3), (f4, ∂tf4),)
  dtf(t) = x -> ForwardDiff.derivative(t -> f(t)(x), t)

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂t(f)(tv)(xv) ≈ ∂tf(tv)(xv)
  @test ∂t(f)(tv)(xv) ≈ dtf(tv)(xv)

  F = TimeSpaceFunction(f)
  @test F(tv)(xv) ≈ f(tv)(xv)
  @test ∂t(F)(tv)(xv) ≈ ∂tf(tv)(xv)
end

# First time derivative, vector-valued
f1(t) = x -> VectorValue(5 * x[1] * x[2], x[2]^2 * t^3)
∂tf1(t) = x -> VectorValue(zero(x[1]), x[2]^2 * 3 * t^2)

f2(t) = x -> VectorValue(x[1]^2, zero(x[2]))
∂tf2(t) = x -> VectorValue(zero(x[1]), zero(x[2]))

f3(t) = x -> VectorValue(x[1]^2, t)
∂tf3(t) = x -> VectorValue(zero(x[1]), one(t))

for (f, ∂tf) in ((f1, ∂tf1), (f2, ∂tf2), (f3, ∂tf3),)
  dtf(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(f(t)(x)), t))

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂t(f)(tv)(xv) ≈ ∂tf(tv)(xv)
  @test ∂t(f)(tv)(xv) ≈ dtf(tv)(xv)

  F = TimeSpaceFunction(f)
  @test F(tv)(xv) ≈ f(tv)(xv)
  @test ∂t(F)(tv)(xv) ≈ ∂tf(tv)(xv)
end

# First time derivative, tensor-valued
f1(t) = x -> TensorValue(x[1] * t, x[2] * t, x[1] * x[2], x[1] * t^2)
∂tf1(t) = x -> TensorValue(x[1], x[2], zero(x[1]), 2 * x[1] * t)

for (f, ∂tf) in ((f1, ∂tf1),)
  dtf(t) = x -> TensorValue(ForwardDiff.derivative(t -> get_array(f(t)(x)), t))

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂t(f)(tv)(xv) ≈ ∂tf(tv)(xv)
  @test ∂t(f)(tv)(xv) ≈ dtf(tv)(xv)

  F = TimeSpaceFunction(f)
  @test F(tv)(xv) ≈ f(tv)(xv)
  @test ∂t(F)(tv)(xv) ≈ ∂tf(tv)(xv)
end

# Spatial derivatives
ft(t) = x -> x[1]^2 * t + x[2]
f = TimeSpaceFunction(ft)

tv = rand(Float64)
xv = Point(rand(Float64, 2)...)
@test ∇(f)(tv)(xv) ≈ Point(2 * xv[1] * tv, 1.0)

# Second time derivative, scalar-valued
f1(t) = x -> t^2
∂tf1(t) = x -> 2 * t
∂ttf1(t) = x -> 2 * one(t)

f2(t) = x -> x[1] * t^2
∂tf2(t) = x -> 2 * x[1] * t
∂ttf2(t) = x -> 2 * x[1]

for (f, ∂tf, ∂ttf) in ((f1, ∂tf1, ∂ttf1), (f2, ∂tf2, ∂ttf2),)
  dtf(t) = x -> ForwardDiff.derivative(t -> f(t)(x), t)
  dttf(t) = x -> ForwardDiff.derivative(t -> dtf(t)(x), t)

  tv = rand(Float64)
  xv = Point(rand(Float64, 2)...)
  @test ∂tt(f)(tv)(xv) ≈ ∂ttf(tv)(xv)
  @test ∂tt(f)(tv)(xv) ≈ dttf(tv)(xv)

  F = TimeSpaceFunction(f)
  @test F(tv)(xv) ≈ f(tv)(xv)
  @test ∂tt(F)(tv)(xv) ≈ ∂ttf(tv)(xv)
end

end # module TimeDerivativesTests
