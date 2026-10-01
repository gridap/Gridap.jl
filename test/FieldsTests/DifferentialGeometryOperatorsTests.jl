module DifferentialGeometryOperatorsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.TensorValues: get_indep_components
using Gridap.Fields

using StaticArrays
using ForwardDiff

const cov2 = Cov(2)

_flat(a::Number)     = SVector(a)
_flat(a::MultiValue) = SVector(Tuple(get_array(a)))
_fullsize(a::Number)     = ()
_fullsize(a::MultiValue) = size(get_array(a))

function ref_partials(f, x0::SVector{D}, idxs::Tuple) where D
  isempty(idxs) && return _flat(f(Point(x0)))
  T = eltype(x0)
  e = SVector{D}(ntuple(k -> k == first(idxs) ? one(T) : zero(T), D))
  ForwardDiff.derivative(t -> ref_partials(f, x0 + t*e, Base.tail(idxs)), zero(T))
end

function test_partial_derivatives(f, x::Point{D}, Nmax) where D
  fx = f(x)
  sz = _fullsize(fx)
  for N in 1:Nmax
    pd = partial_derivatives(f,x,Val(N))
    @test typeof(pd) === derivatives_type(typeof(fx),x,Val(N))
    @test size(pd) === (ntuple(_->Cov(D),N)..., size(fx)...)
    @test (@inferred partial_derivatives(f,x,Val(N))) === pd
    for I in CartesianIndices(ntuple(_->D,N))
      r = ref_partials(f,get_array(x),Tuple(I))
      for (n,J) in enumerate(CartesianIndices(sz))
        @test pd[Tuple(I)...,Tuple(J)...] ≈ r[n] atol=1e-9
      end
    end
  end
end

x2 = Point(0.3,-0.7)
x3 = Point(0.3,-0.7,0.2)



# scalar valued, cross checked against ForwardDiff's own gradient and hessian
fs(x) = sin(x[1])*exp(x[2]) + x[1]^3*x[2]

test_partial_derivatives(fs, x2, 2)
@test get_array(partial_derivatives(fs,x2,Val(1))) ≈ ForwardDiff.gradient(fs∘Point,get_array(x2))
@test get_array(partial_derivatives(fs,x2,Val(2))) ≈ ForwardDiff.hessian(fs∘Point,get_array(x2))

fi(x) = 3

@test partial_derivatives(fs,x2,Val(0)) === fs(x2)
@test partial_derivatives(fi,x2,Val(0)) === 3.0
@test typeof(partial_derivatives(fi,x2,Val(0))) === derivatives_type(Int,x2,Val(0))

# A vector valued function has the index layout of `gradient`, ∂u[i,j] = ∂u_j/∂xⁱ
fv(x) = VectorValue(sin(x[1])*x[2], exp(x[2])+x[1]^3)

@test get_array(partial_derivatives(fv,x2,Val(1))) ≈ get_array(gradient(fv,x2))
test_partial_derivatives(fv, x2, 2)

# Every tensor family, including those whose stored components are dependent
w(l,x) = sin((l+1)*x[1]) + l*x[2]^2 + x[1]*x[2]/(l+1)
mk(V,x) = V(ntuple(l -> w(l,x), Val(num_indep_components(V))))

for V in (VectorValue{2}, CoVectorValue{cov2}, TensorValue{2,2},
          SymTensorValue{2}, SymTracelessTensorValue{2},
          SkewSymTensorValue{2}, ExteriorFormValue{1,2},
          ExteriorFormValue{2,2}, ContraVectorValue{Con(2)})
  test_partial_derivatives(x -> mk(V,x), x2, 2)
end

for V in (VectorValue{3}, SymTensorValue{3}, ExteriorFormValue{2,3},
          ThirdOrderTensorValue{3,3,3})
  test_partial_derivatives(x -> mk(V,x), x3, 2)
end

# A value type without independent component has constant, zero components
f0(x) = zero(ExteriorFormValue{3,2,Float64,0})

@test iszero(partial_derivatives(f0,x2,Val(1)))
@test typeof(partial_derivatives(f0,x2,Val(1))) ===
  derivatives_type(ExteriorFormValue{3,2,Float64,0},x2,Val(1))

@test partial_derivatives(fv,Val(2))(x2) === partial_derivatives(fv,x2,Val(2))
@test partial_derivatives(fv)(x2) === partial_derivatives(fv,x2,Val(1))

# Partials of a tensor that is not an alternating form carry neither exterior
# derivative nor codifferential, and the Lie derivative is declared, not implemented
pdv = partial_derivatives(fv,x2,Val(1))

@test_throws Exception partials_to_exteriorderivative(pdv)
@test_throws Exception partials_to_liederivative(fv(x2),fv(x2),pdv,pdv)
@test_throws Exception partials_to_codifferential(pdv)

# Exterior derivative of a K-form from its first order partial derivatives

d(f) = y -> partials_to_exteriorderivative(partial_derivatives(f,y,Val(1)))

function test_exterior_derivative(K, x::Point{D}) where D
  V = ExteriorFormValue{K,Cov(D)}
  f = y -> mk(V,y)
  ∂ω = partial_derivatives(f,x,Val(1))
  dω = partials_to_exteriorderivative(∂ω)

  @test typeof(dω) === ExteriorFormValue{K+1,Cov(D),Float64,binomial(D,K+1)}
  @test (@inferred partials_to_exteriorderivative(∂ω)) === dω
  # dω is the alternation of ∂ω scaled by K+1, which reads ∂ω by full indices
  @test dω ≈ (K+1)*alternation(∂ω)
  # the mixed partials are symmetric, so d∘d vanishes exactly
  @test iszero(d(d(f))(x))
end

for K in 0:2
  test_exterior_derivative(K, x2)
end
for K in 0:3
  test_exterior_derivative(K, x3)
end

# A 1-form has the exterior derivative of the vector field proxying it
fv3(y) = VectorValue(sin(y[1])*y[2], exp(y[2])+y[1]^3, y[3]*y[1]^2+y[2])
proxied_1form(fvec, ::Val{D}) where D = y -> ExteriorFormValue{1,Cov(D)}(Tuple(fvec(y)))

@test d(proxied_1form(fv, Val(2)))(x2) ≈ grad_to_2form(gradient(fv,x2))
@test d(proxied_1form(fv3,Val(3)))(x3) ≈ grad_to_2form(gradient(fv3,x3))

# A form with no component differentiates to one with no component
fe(y) = zero(ExteriorFormValue{3,Cov(2),Float64,0})
@test iszero(d(fe)(x2))
@test typeof(d(fe)(x2)) === ExteriorFormValue{4,Cov(2),Float64,0}

# The exterior derivative is defined on forms, not on K-vectors
pdc = partial_derivatives(y -> mk(ExteriorVectorValue{1,Con(2)},y), x2, Val(1))
@test_throws Exception partials_to_exteriorderivative(pdc)
@test_throws ArgumentError partials_to_exteriorderivative(
  SymNTensorValue{Tuple{Con(2)}}(1.0,2.0))

# A scalar function is a 0-form proxied by a Number: its partials are a bare
# order 1 SymNTensorValue, and d of it is the gradient as a 1-form
for (f, x) in ((fs, x2), (y -> sin(y[1])*exp(y[2]) + y[1]^3*y[2]*y[3], x3))
  D = length(x)
  ∂f = partial_derivatives(f,x,Val(1))
  @test ∂f isa SymNTensorValue{Tuple{Cov(D)},Float64,1,D}
  df = partials_to_exteriorderivative(∂f)
  @test typeof(df) === ExteriorFormValue{1,Cov(D),Float64,D}
  @test (@inferred partials_to_exteriorderivative(∂f)) === df
  @test Tuple(df) == Tuple(gradient(f,x))
  # d∘d = 0 on the 1-form it produces
  @test iszero(d(y -> partials_to_exteriorderivative(partial_derivatives(f,y,Val(1))))(x))
end

# Codifferential of a K-form from its first order partial derivatives

δ(f) = y -> partials_to_codifferential(partial_derivatives(f,y,Val(1)))

function test_codifferential(K, x::Point{D}) where D
  V = ExteriorFormValue{K,Cov(D)}
  f = y -> mk(V,y)
  ∂ω = partial_derivatives(f,x,Val(1))
  δω = partials_to_codifferential(∂ω)

  @test typeof(δω) === ExteriorFormValue{K-1,Cov(D),Float64,binomial(D,K-1)}
  @test (@inferred partials_to_codifferential(∂ω)) === δω
  # the derivation index contracts with the first index of ω, reading ∂ω by full indices
  for (n,I) in enumerate(sorted_combinations(D,K-1))
    @test indep_comp_getindex(δω,n) ≈ -sum(j -> ∂ω[j,j,I...], 1:D)
  end
  # δ is ±⋆d⋆ in Euclidean space
  @test δω ≈ (isodd(D*(K-1)+1) ? -1 : 1) * ⋆ᴱ(d(y -> ⋆ᴱ(f(y)))(x))
  # the mixed partials are symmetric, so δ∘δ vanishes exactly where it is defined
  K ≥ 2 && @test iszero(δ(δ(f))(x))
end

for K in 1:2
  test_codifferential(K, x2)
end
for K in 1:3
  test_codifferential(K, x3)
end

# The codifferential of a 0-form would be a (-1)-form
f0form(y) = mk(ExteriorFormValue{0,Cov(2)},y)
@test_throws Exception partials_to_codifferential(partial_derivatives(f0form,x2,Val(1)))

# A form with no component has a vanishing codifferential, its sole component
# summing no partial derivative at all
@test iszero(δ(fe)(x2))
@test typeof(δ(fe)(x2)) === ExteriorFormValue{2,Cov(2),Float64,1}

f4e(y) = zero(ExteriorFormValue{4,Cov(2),Float64,0})
@test typeof(δ(f4e)(x2)) === ExteriorFormValue{3,Cov(2),Float64,0}

# The derivation indices may be held by any order 1 tensor, not only by the
# SymNTensorValue that partial_derivatives produces
∂2form = ExteriorFormValue{2,Cov(3)}(5.,6.,7.)
∂sym   = product_tensor(SymNTensorValue{Tuple{Cov(3)}}(1.,2.,3.), ∂2form)
∂covec = product_tensor(CoVectorValue{Cov(3)}(1.,2.,3.), ∂2form)
@test partials_to_exteriorderivative(∂covec) === partials_to_exteriorderivative(∂sym)
@test partials_to_codifferential(∂covec) === partials_to_codifferential(∂sym)

# The 0-form method takes any order 1 covariant tensor as well
@test partials_to_exteriorderivative(CoVectorValue{Cov(2)}(1.,2.)) ===
  partials_to_exteriorderivative(SymNTensorValue{Tuple{Cov(2)}}(1.,2.))
@test_throws ArgumentError partials_to_exteriorderivative(VectorValue(1.,2.))

# The codifferential is defined on forms, not on K-vectors
@test_throws Exception partials_to_codifferential(pdc)

# Lie derivative along a contravariant vector field

function test_liederivative_form(K, x::Point{D}) where D
  V = ExteriorFormValue{K,Cov(D)}
  uF = y -> ContraVectorValue{Con(D)}(ntuple(l -> w(l,y) + y[1], D))
  ωF = y -> mk(V,y)
  u, ω = uF(x), ωF(x)
  ∂u, ∂ω = partial_derivatives(uF,x,Val(1)), partial_derivatives(ωF,x,Val(1))
  Lω = partials_to_liederivative(u, ω, ∂u, ∂ω)

  @test typeof(Lω) === ExteriorFormValue{K,Cov(D),Float64,binomial(D,K)}
  @test (@inferred partials_to_liederivative(u, ω, ∂u, ∂ω)) === Lω
  # Cartan's magic formula 𝓛ᵤω = d(ι_u ω) + ι_u(dω), undefined for a 0-form
  if K ≥ 1
    d_iuω = partials_to_exteriorderivative(
      partial_derivatives(y -> interior_product(uF(y), ωF(y)), x, Val(1)))
    @test Lω ≈ d_iuω + interior_product(u, partials_to_exteriorderivative(∂ω))
  end
end

for K in 0:2
  test_liederivative_form(K, x2)
end
for K in 0:3
  test_liederivative_form(K, x3)
end

uvF(y) = ContraVectorValue{Con(2)}(sin(y[1])+y[2]^2, exp(y[2])*y[1])
vvF(y) = ContraVectorValue{Con(2)}(cos(y[2])*y[1], y[1]^3+y[2])
αF(y)  = CoVectorValue{Cov(2)}(y[1]*y[2]+sin(y[2]), cos(y[1])+y[2]^2)
∂uv, ∂vv, ∂α = partial_derivatives(uvF,x2,Val(1)), partial_derivatives(vvF,x2,Val(1)),
                partial_derivatives(αF,x2,Val(1))
Luv = partials_to_liederivative(uvF(x2), vvF(x2), ∂uv, ∂vv)
Lvu = partials_to_liederivative(vvF(x2), uvF(x2), ∂vv, ∂uv)
Luα = partials_to_liederivative(uvF(x2), αF(x2), ∂uv, ∂α)

@test typeof(Luv) === ContraVectorValue{Con(2),Float64,2}
@test typeof(Luα) === CoVectorValue{Cov(2),Float64,2}
# 𝓛ᵤv is the Lie bracket [u,v], which is antisymmetric
@test all(Tuple(Luv) .≈ .-Tuple(Lvu))
# ⟨α,v⟩ obeys Leibniz along u, which fixes the sign of both Lie derivatives
∂pair = partial_derivatives(y -> interior_product(vvF(y), αF(y)), x2, Val(1))
@test sum(Tuple(∂pair) .* Tuple(uvF(x2))) ≈
  interior_product(vvF(x2), Luα)[1] + interior_product(Luv, αF(x2))[1]

# An order 1 covariant tensor that is not a 1-form takes the same covector formula
βF(y) = SymNTensorValue{Tuple{Cov(2)}}(y[1]*y[2]+sin(y[2]), cos(y[1])+y[2]^2)
Luβ = partials_to_liederivative(uvF(x2), βF(x2), ∂uv, partial_derivatives(βF,x2,Val(1)))
@test all(Tuple(Luβ) .≈ Tuple(Luα))

# There is nothing to flow along without a contravariant first argument
@test_throws Exception partials_to_liederivative(αF(x2), vvF(x2), ∂α, ∂vv)
@test_throws Exception partials_to_liederivative(αF(x2), αF(x2), ∂α, ∂α)

# Lie derivative of a scalar function, the derivative of f along the flow of u
function test_liederivative_scalar(f, x::Point{D}) where D
  uF = y -> ContraVectorValue{Con(D)}(ntuple(l -> w(l,y) + y[1], D))
  u, ∂u = uF(x), partial_derivatives(uF,x,Val(1))
  ∂f = partial_derivatives(f,x,Val(1))
  Lf = partials_to_liederivative(u, f(x), ∂u, ∂f)

  @test Lf isa Float64
  @test (@inferred partials_to_liederivative(u, f(x), ∂u, ∂f)) === Lf
  # 𝓛ᵤf = uʲ ∂ⱼf is the directional derivative of f along u
  @test Lf ≈ sum(Tuple(u) .* Tuple(gradient(f,x)))
  # a scalar function is a 0-form proxied by a Number, both give the same value
  f0F = y -> ExteriorFormValue{0,Cov(D)}(f(y))
  Lf0 = partials_to_liederivative(u, f0F(x), ∂u, partial_derivatives(f0F,x,Val(1)))
  @test typeof(Lf0) === ExteriorFormValue{0,Cov(D),Float64,1}
  @test Lf0[1] ≈ Lf
end

test_liederivative_scalar(fs, x2)
test_liederivative_scalar(y -> sin(y[1])*exp(y[2]) + y[1]^3*y[2]*y[3], x3)

# There is nothing to flow along without a contravariant first argument, and the
# derivation index of the partial derivatives must be covariant
∂fs = partial_derivatives(fs,x2,Val(1))
@test_throws ArgumentError partials_to_liederivative(
  αF(x2), fs(x2), TensorValue{Cov(2),Cov(2)}(1.,2.,3.,4.), ∂fs)
@test_throws ArgumentError partials_to_liederivative(
  uvF(x2), fs(x2), TensorValue{Con(2),Con(2)}(1.,2.,3.,4.),
  SymNTensorValue{Tuple{Con(2)}}(1.,2.))

# A `MultiValue` is a `Number`, so the scalar method must stay more specific than
# the not-implemented fallback, which keeps catching the shapes it does not cover
@test_throws ErrorException partials_to_liederivative(
  uvF(x2), TensorValue{2,2}(1.,2.,3.,4.), ∂uv, ∂fs)

end # module
