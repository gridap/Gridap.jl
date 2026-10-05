module HigherOrderGradientsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: seed_point, taylor_tensor, nth_gradient_type, push_∇ⁿ

# A scalar polynomial field p(x) = Σ_α c_α x^α with no hand-written gradient: all its
# derivatives come from the AD default. Its derivatives in closed form are the reference.
struct PolyField{D} <: Field
  coeffs::Vector{Pair{NTuple{D,Int},Float64}}
end

function Fields.evaluate!(cache, f::PolyField{D}, x::Point) where D
  sum(c * prod(x[d]^α[d] for d in 1:D) for (α, c) in f.coeffs)
end

function dmonomial(α::NTuple{D,Int}, dirs, x) where D
  β = collect(α)
  c = 1.0
  for d in dirs
    β[d] == 0 && return 0.0
    c *= β[d]
    β[d] -= 1
  end
  c * prod(x[d]^β[d] for d in 1:D)
end

function exact_gradient(f::PolyField{D}, x::Point{D}, N) where D
  G = nth_gradient_type(Float64, x, Val(N))
  vals = [sum(c * dmonomial(α, Tuple(I), x) for (α, c) in f.coeffs) for I in CartesianIndices(ntuple(_ -> D, N))]
  G(Tuple(vals))
end

maxdiff(a::MultiValue, b::MultiValue) = maximum(abs, Tuple(a - b))
maxdiff(a::AbstractArray, b::AbstractArray) = maximum(maxdiff.(a, b))

# Gradient types

x2 = Point(0.3, 0.7)
@test nth_gradient_type(Float64, x2, Val(0)) == Float64
@test nth_gradient_type(Float64, x2, Val(1)) == VectorValue{2,Float64}
@test nth_gradient_type(Float64, x2, Val(2)) == TensorValue{2,2,Float64,4}
@test nth_gradient_type(Float64, x2, Val(3)) == ThirdOrderTensorValue{2,2,2,Float64,8}
@test nth_gradient_type(Float64, x2, Val(4)) == HighOrderTensorValue{Tuple{2,2,2,2},Float64,4,16}
@test nth_gradient_type(VectorValue{2,Float64}, x2, Val(2)) == ThirdOrderTensorValue{2,2,2,Float64,8}

# Seeding and extraction

f(x) = x[1]^3 * x[2]^2
v = f(seed_point(x2, Val(2)))
H = taylor_tensor(TensorValue{2,2,Float64,4}, v, Val(2), Val(2))
@test H ≈ TensorValue(6x2[1]*x2[2]^2, 6x2[1]^2*x2[2], 6x2[1]^2*x2[2], 2x2[1]^3)

# Single field: FieldGradient{N} by AD, any N

p2 = PolyField{2}([(3,3)=>1.0, (2,1)=>-2.0, (1,3)=>0.5, (4,0)=>1.5, (0,2)=>1.0])
for N in 1:5
  @test maxdiff(evaluate(gradient(p2, Val(N)), x2), exact_gradient(p2, x2, N)) < 1e-10
end
@test evaluate(∇(p2), x2) ≈ exact_gradient(p2, x2, 1)
@test evaluate(∇∇(p2), x2) ≈ exact_gradient(p2, x2, 2)
@test gradient(p2, Val(3)) isa FieldGradient{3}
@test gradient(gradient(p2, Val(2)), Val(2)) isa FieldGradient{4}

x3 = Point(0.3, 0.7, 0.2)
p3 = PolyField{3}([(2,2,2)=>1.0, (2,1,0)=>-1.0, (0,1,3)=>0.5])
for N in 1:3
  @test maxdiff(evaluate(gradient(p3, Val(N)), x3), exact_gradient(p3, x3, N)) < 1e-10
end

# Arrays of fields: FieldGradientArray{N} by AD, the whole array evaluated at once

fields = [p2, PolyField{2}([(1,1)=>2.0]), PolyField{2}([(0,4)=>-1.0, (2,0)=>3.0])]
xs = [Point(0.1, 0.2), Point(0.5, 0.9), Point(0.3, 0.7), Point(1.0, -0.5)]
for N in 1:4
  g = Broadcasting(NthGradient{N}())(fields)
  @test g isa FieldGradientArray{N}
  r = evaluate(g, xs)
  @test size(r) == (length(xs), length(fields))
  @test maxdiff(r, [exact_gradient(fi, xi, N) for xi in xs, fi in fields]) < 1e-10
  test_map(r, g, xs; cmp = (a, b) -> maxdiff(a, b) < 1e-12)
end
@test Broadcasting(∇)(Broadcasting(NthGradient{3}())(fields)) isa FieldGradientArray{4}
@test Broadcasting(NthGradient{2}())(Broadcasting(NthGradient{3}())(fields)) isa FieldGradientArray{5}

# Linear combinations of fields differentiate the underlying fields

values = [1.0 2.0; -1.0 0.5; 0.0 3.0]
lc = linear_combination(values, fields)
r = evaluate(Broadcasting(NthGradient{3}())(lc), xs)
r_exact = [sum(values[i,j] * exact_gradient(fields[i], xk, 3) for i in 1:3) for xk in xs, j in 1:2]
@test maxdiff(r, r_exact) < 1e-10

# Pushforward by an affine map: push_∇ and push_∇∇ at N = 1, 2

Jt_inv = TensorValue(0.5, 0.1, -0.2, 2.0)
g1 = VectorValue(1.0, -3.0)
@test push_∇ⁿ(g1, Jt_inv, Val(1)) ≈ Jt_inv ⋅ g1
g2 = TensorValue(1.0, 2.0, 2.0, -1.0)
@test push_∇ⁿ(g2, Jt_inv, Val(2)) ≈ push_∇∇(g2, Jt_inv)

end # module
