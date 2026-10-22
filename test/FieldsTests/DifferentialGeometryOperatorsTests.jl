module DifferentialGeometryOperatorsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.TensorValues: _sqrt_abs_det

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
    pd = PartialDerivatives{N}()(f,x)
    @test typeof(pd) === derivatives_type(typeof(fx),x,Val(N))
    @test size(pd)[1:N] === ntuple(_->Cov(D),N)
    @test (@inferred PartialDerivatives{N}()(f,x)) === pd
    for I in CartesianIndices(ntuple(_->D,N))
      r = ref_partials(f,get_array(x),Tuple(I))
      for (n,J) in enumerate(CartesianIndices(sz))
        @test pd[Tuple(I)...,Tuple(J)...] ≈ r[n] atol=1e-9
      end
    end
  end
end

# ── PartialDerivatives ────────────────────────────────────────────────────────
# The named orders, and the constructor guard on the order itself

@test partial_derivatives        === PartialDerivatives{1}()
@test second_partial_derivatives === PartialDerivatives{2}()
@test_throws ArgumentError PartialDerivatives{-1}()
@test_throws ArgumentError PartialDerivatives{1.0}()
@test_throws ArgumentError PartialDerivatives{(1,2)}()

x2 = Point(0.3,-0.7)
x3 = Point(0.3,-0.7,0.2)

# scalar valued, cross checked against ForwardDiff's own gradient and hessian
fs(x) = sin(x[1])*exp(x[2]) + x[1]^3*x[2]

test_partial_derivatives(fs, x2, 2)
@test get_array(partial_derivatives(fs,x2)) ≈ ForwardDiff.gradient(fs∘Point,get_array(x2))
@test get_array(second_partial_derivatives(fs,x2)) ≈ ForwardDiff.hessian(fs∘Point,get_array(x2))

fi(x) = 3

@test PartialDerivatives{0}()(fs,x2) === fs(x2)
@test PartialDerivatives{0}()(fi,x2) === 3.0
@test typeof(PartialDerivatives{0}()(fi,x2)) === derivatives_type(Int,x2,Val(0))

# A vector valued function has the index layout of `gradient`, ∂u[i,j] = ∂u_j/∂xⁱ
fv(x) = VectorValue(sin(x[1])*x[2], exp(x[2])+x[1]^3)

@test get_array(partial_derivatives(fv,x2)) ≈ get_array(gradient(fv,x2))
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

@test iszero(partial_derivatives(f0,x2))
@test typeof(partial_derivatives(f0,x2)) ===
  derivatives_type(ExteriorFormValue{3,2,Float64,0},x2,Val(1))

@test second_partial_derivatives(fv)(x2) === second_partial_derivatives(fv,x2)
@test partial_derivatives(fv)(x2) === partial_derivatives(fv,x2)

# Partials of a tensor that is not an alternating form carry neither exterior
# derivative nor codifferential, and the Lie derivative is declared, not implemented
pdv = partial_derivatives(fv,x2)

@test_throws Exception partials_to_exteriorderivative(pdv)
@test_throws Exception partials_to_liederivative(fv(x2),fv(x2),pdv,pdv)
@test_throws Exception partials_to_flat_codifferential(pdv)

# Exterior derivative of a K-form from its first order partial derivatives

d(f) = y -> partials_to_exteriorderivative(partial_derivatives(f,y))

function test_exterior_derivative(K, x::Point{D}) where D
  V = ExteriorFormValue{K,Cov(D)}
  f = y -> mk(V,y)
  ∂ω = partial_derivatives(f,x)
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
pdc = partial_derivatives(y -> mk(ExteriorVectorValue{1,Con(2)},y), x2)
@test_throws Exception partials_to_exteriorderivative(pdc)
@test_throws ArgumentError partials_to_exteriorderivative(
  SymNTensorValue{Tuple{Con(2)}}(1.0,2.0))

# A scalar function is a 0-form proxied by a Number: its partials are an order 1
# SymNTensorValue tensored with the order 0 ScalarValue carrying the value, so
# they are an order 1 covariant tensor, and d of it is the gradient as a 1-form
for (f, x) in ((fs, x2), (y -> sin(y[1])*exp(y[2]) + y[1]^3*y[2]*y[3], x3))
  D = length(x)
  ∂f = partial_derivatives(f,x)
  @test ∂f isa ProductTensorValue{SymNTensorValue{Tuple{Cov(D)}},ScalarValue,Float64,Tuple{Cov(D)},1,D}
  @test ∂f isa MultiValue{Tuple{Cov(D)},Float64,1,D}
  df = partials_to_exteriorderivative(∂f)
  @test typeof(df) === ExteriorFormValue{1,Cov(D),Float64,D}
  @test (@inferred partials_to_exteriorderivative(∂f)) === df
  @test Tuple(df) == Tuple(gradient(f,x))
  # d∘d = 0 on the 1-form it produces
  @test iszero(d(y -> partials_to_exteriorderivative(partial_derivatives(f,y)))(x))
end

# Codifferential of a K-form from its first order partial derivatives

δ(f) = y -> partials_to_flat_codifferential(partial_derivatives(f,y))

function test_codifferential(K, x::Point{D}) where D
  V = ExteriorFormValue{K,Cov(D)}
  f = y -> mk(V,y)
  ∂ω = partial_derivatives(f,x)
  δω = partials_to_flat_codifferential(∂ω)

  @test typeof(δω) === ExteriorFormValue{K-1,Cov(D),Float64,binomial(D,K-1)}
  @test (@inferred partials_to_flat_codifferential(∂ω)) === δω
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
@test_throws Exception partials_to_flat_codifferential(partial_derivatives(f0form,x2))

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
@test partials_to_flat_codifferential(∂covec) === partials_to_flat_codifferential(∂sym)

# The 0-form method takes any order 1 covariant tensor as well
@test partials_to_exteriorderivative(CoVectorValue{Cov(2)}(1.,2.)) ===
  partials_to_exteriorderivative(SymNTensorValue{Tuple{Cov(2)}}(1.,2.))
@test_throws ArgumentError partials_to_exteriorderivative(VectorValue(1.,2.))

# The codifferential is defined on forms, not on K-vectors
@test_throws Exception partials_to_flat_codifferential(pdc)

# Lie derivative along a contravariant vector field

function test_liederivative_form(K, x::Point{D}) where D
  uF = y -> ContraVectorValue{Con(D)}(ntuple(l -> w(l,y) + y[1], D))
  ωF = y -> mk(ExteriorFormValue{K,Cov(D)}, y)
  u, ω = uF(x), ωF(x)
  ∂u, ∂ω = partial_derivatives(uF,x), partial_derivatives(ωF,x)
  Lω = partials_to_liederivative(u, ω, ∂u, ∂ω)

  @test typeof(Lω) === ExteriorFormValue{K,Cov(D),Float64,binomial(D,K)}
  @test (@inferred partials_to_liederivative(u, ω, ∂u, ∂ω)) === Lω
  # Cartan's magic formula 𝓛ᵤω = d(ι_u ω) + ι_u(dω), undefined for a 0-form
  if K ≥ 1
    d_iuω = partials_to_exteriorderivative(
      partial_derivatives(y -> interior_product(uF(y), ωF(y)), x))
    @test Lω ≈ d_iuω + interior_product(u, partials_to_exteriorderivative(∂ω))
  end
end

for K in 0:2
  test_liederivative_form(K, x2)
end
for K in 0:3
  test_liederivative_form(K, x3)
end

# K-vector
function test_liederivative_vector(K, x::Point{D}) where D
  uF = y -> ContraVectorValue{Con(D)}(ntuple(l -> w(l,y) + y[1], D))
  VF = y -> mk(ExteriorVectorValue{K,Con(D)}, y)
  u, ∂u = uF(x), partial_derivatives(uF,x)
  ∂V = partial_derivatives(VF,x)
  LV = partials_to_liederivative(u, VF(x), ∂u, ∂V)

  @test typeof(LV) === ExteriorVectorValue{K,Con(D),Float64,binomial(D,K)}
  @test (@inferred partials_to_liederivative(u, VF(x), ∂u, ∂V)) === LV
  # A K = 0 pairing is a product of scalars, which carries no index to transport
  if K ≥ 1
    ωF = y -> mk(ExteriorFormValue{K,Cov(D)}, y)
    Lω = partials_to_liederivative(u, ωF(x), ∂u, partial_derivatives(ωF,x))
    pairing(a,b) = contracted_product(Val(K), a, b)
    ∂pairing = partial_derivatives(y -> pairing(ωF(y), VF(y)), x)
    @test sum(Tuple(u) .* Tuple(∂pairing)) ≈ pairing(Lω, VF(x)) + pairing(ωF(x), LV)
  end
end

for K in 0:2
  test_liederivative_vector(K, x2)
end
for K in 0:3
  test_liederivative_vector(K, x3)
end

uvF(y) = ContraVectorValue{Con(2)}(sin(y[1])+y[2]^2, exp(y[2])*y[1])
vvF(y) = ContraVectorValue{Con(2)}(cos(y[2])*y[1], y[1]^3+y[2])
αF(y)  = CoVectorValue{Cov(2)}(y[1]*y[2]+sin(y[2]), cos(y[1])+y[2]^2)
∂uv, ∂vv, ∂α = partial_derivatives(uvF,x2), partial_derivatives(vvF,x2),
                partial_derivatives(αF,x2)
Luv = partials_to_liederivative(uvF(x2), vvF(x2), ∂uv, ∂vv)
Lvu = partials_to_liederivative(vvF(x2), uvF(x2), ∂vv, ∂uv)
Luα = partials_to_liederivative(uvF(x2), αF(x2), ∂uv, ∂α)

@test typeof(Luv) === ContraVectorValue{Con(2),Float64,2}
@test typeof(Luα) === CoVectorValue{Cov(2),Float64,2}
# 𝓛ᵤv is the Lie bracket [u,v], which is antisymmetric
@test all(Tuple(Luv) .≈ .-Tuple(Lvu))
# ⟨α,v⟩ obeys Leibniz along u, which fixes the sign of both Lie derivatives
∂pair = partial_derivatives(y -> interior_product(vvF(y), αF(y)), x2)
@test sum(Tuple(∂pair) .* Tuple(uvF(x2))) ≈
  interior_product(vvF(x2), Luα)[1] + interior_product(Luv, αF(x2))[1]

# An order 1 covariant tensor that is not a 1-form takes the same covector formula
βF(y) = SymNTensorValue{Tuple{Cov(2)}}(y[1]*y[2]+sin(y[2]), cos(y[1])+y[2]^2)
Luβ = partials_to_liederivative(uvF(x2), βF(x2), ∂uv, partial_derivatives(βF,x2))
@test all(Tuple(Luβ) .≈ Tuple(Luα))

# There is nothing to flow along without a contravariant first argument
@test_throws Exception partials_to_liederivative(αF(x2), vvF(x2), ∂α, ∂vv)
@test_throws Exception partials_to_liederivative(αF(x2), αF(x2), ∂α, ∂α)

# Lie derivative of a scalar function, the derivative of f along the flow of u
function test_liederivative_scalar(f, x::Point{D}) where D
  uF = y -> ContraVectorValue{Con(D)}(ntuple(l -> w(l,y) + y[1], D))
  u, ∂u = uF(x), partial_derivatives(uF,x)
  ∂f = partial_derivatives(f,x)
  Lf = partials_to_liederivative(u, f(x), ∂u, ∂f)

  @test Lf isa Float64
  @test (@inferred partials_to_liederivative(u, f(x), ∂u, ∂f)) === Lf
  # 𝓛ᵤf = uʲ ∂ⱼf is the directional derivative of f along u
  @test Lf ≈ sum(Tuple(u) .* Tuple(gradient(f,x)))
  # a scalar function is a 0-form proxied by a Number, both give the same value
  f0F = y -> ExteriorFormValue{0,Cov(D)}(f(y))
  Lf0 = partials_to_liederivative(u, f0F(x), ∂u, partial_derivatives(f0F,x))
  @test typeof(Lf0) === ExteriorFormValue{0,Cov(D),Float64,1}
  @test Lf0[1] ≈ Lf
end

test_liederivative_scalar(fs, x2)
test_liederivative_scalar(y -> sin(y[1])*exp(y[2]) + y[1]^3*y[2]*y[3], x3)

# A tensor with no (or unknown) index symmetry takes the general formula, each of
# its indices transported according to its own variance.

# 𝓛ᵤδ = 0 for the identity (1,1) tensor: its covariant index is transported
# exactly against its contravariant one, which fixes the sign of both
δF(y) = TensorValue{Con(2),Cov(2)}(1.,0.,0.,1.)
∂δ = partial_derivatives(δF,x2)
Lδ = partials_to_liederivative(uvF(x2), δF(x2), ∂uv, ∂δ)
@test typeof(Lδ) === TensorValue{Con(2),Cov(2),Float64,4}
@test (@inferred partials_to_liederivative(uvF(x2), δF(x2), ∂uv, ∂δ)) === Lδ
@test all(iszero, Tuple(Lδ))

# On the flat metric 𝓛ᵤg is the Killing operator ∂u + ∂uᵀ, so it vanishes for a
# rotation field and keeps the symmetry of `g`
gF(y) = SymTensorValue{cov2}(1.,0.,1.)
∂g = partial_derivatives(gF,x2)
Lg = partials_to_liederivative(uvF(x2), gF(x2), ∂uv, ∂g)
@test typeof(Lg) === SymTensorValue{cov2,Float64,3}
@test all(Tuple(Lg) .≈ (2∂uv[1,1], ∂uv[1,2]+∂uv[2,1], 2∂uv[2,2]))
rotF(y) = ContraVectorValue{Con(2)}(-y[2], y[1])
@test all(iszero, Tuple(partials_to_liederivative(
  rotF(x2), gF(x2), partial_derivatives(rotF,x2), ∂g)))

# There is nothing to flow along without a contravariant first argument, and the
# derivation index of the partial derivatives must be covariant
∂fs = partial_derivatives(fs,x2)
@test_throws ArgumentError partials_to_liederivative(
  αF(x2), fs(x2), TensorValue{Cov(2),Cov(2)}(1.,2.,3.,4.), ∂fs)
@test_throws ArgumentError partials_to_liederivative(
  uvF(x2), fs(x2), TensorValue{Con(2),Con(2)}(1.,2.,3.,4.),
  SymNTensorValue{Tuple{Con(2)}}(1.,2.))

# Fallback if unvarianced argument was given.
@test_throws ErrorException partials_to_liederivative(
  uvF(x2), TensorValue{Con(2),Cov(2)}(1.,2.,3.,4.), ∂uv,
  NTensorValue{Tuple{Cov(2),Con(2),Cov(2)}}(1.,2.,3.,4.,5.,6.,7.,8.))

# Christoffel symbols of both kinds, from the partial derivatives of a metric

using LinearAlgebra: I
using Gridap.TensorValues: ProductTensorValue, _indep_comp_cartesian_indices,
  _2d_sym_tensor_index_pair, _2d_sym_tensor_linear_index

# The pair map inverts the index map a `SymTensorValue` is read with, and gives the
# same pairs as the general component machinery, each up to its transposition
for D in 0:6
  for i in 1:D, j in i:D
    @test _2d_sym_tensor_index_pair(D, _2d_sym_tensor_linear_index(D,i,j)) == (i,j)
  end
  D == 0 && continue
  @test [_2d_sym_tensor_index_pair(D,l) for l in 1:(D*(D+1)÷2)] ==
    [minmax(p...) for p in _indep_comp_cartesian_indices(SymTensorValue{Cov(D),Float64})]
end

# Γ_{kij} = (∂ᵢg_{jk} + ∂ⱼg_{ik} - ∂ₖg_{ij})/2, and Γᵏᵢⱼ = g⁻¹ᵏᵐ Γ_{mij}
ref_Christoffel_1(∂g,k,i,j) = (∂g[i,j,k] + ∂g[j,i,k] - ∂g[k,i,j])/2
ref_Christoffel_2(∂g,g_inv,D,k,i,j) =
  sum(g_inv[k,m]*ref_Christoffel_1(∂g,m,i,j) for m in 1:D)

function test_Christoffel(gF, x::Point{D}) where D
  g, ∂g = gF(x), partial_derivatives(gF,x)
  g_inv = inv(g)
  Γ1 = partials_to_Christoffel_1(∂g)
  Γ2 = partials_to_Christoffel_2(∂g, g_inv)

  # the first kind is fully covariant, the second has its first index raised
  @test size(Γ1) === (Cov(D), Cov(D), Cov(D))
  @test size(Γ2) === (Con(D), Cov(D), Cov(D))
  @test num_indep_components(typeof(Γ1)) == D*(D*(D+1)÷2)
  @test (@inferred partials_to_Christoffel_1(∂g)) === Γ1
  @test (@inferred partials_to_Christoffel_2(∂g, g_inv)) === Γ2

  for k in 1:D, i in 1:D, j in 1:D
    @test Γ1[k,i,j] ≈ ref_Christoffel_1(∂g,k,i,j) atol=1e-12
    @test Γ2[k,i,j] ≈ ref_Christoffel_2(∂g,g_inv,D,k,i,j) atol=1e-12
    # both kinds are symmetric in their last two indices
    @test Γ1[k,i,j] === Γ1[k,j,i]
    @test Γ2[k,i,j] === Γ2[k,j,i]
    # lowering the first index of the second kind gives the first kind back
    @test Γ1[k,i,j] ≈ sum(g[k,m]*Γ2[m,i,j] for m in 1:D) atol=1e-12
  end

  # the inverse metric is taken from the metric data in whichever form it comes
  @test partials_to_Christoffel_2(∂g, g) ≈ Γ2
  @test partials_to_Christoffel_2(∂g, parse_metric_data(g_inv)) === Γ2
  @test partials_to_Christoffel_2(∂g, compute_metric_data(g)) ≈ Γ2
  @test (@inferred partials_to_Christoffel_2(∂g, g)) ≈ Γ2

  # the identity metric only raises the first index, leaving the components alone
  ΓI = partials_to_Christoffel_2(∂g, I)
  @test Tuple(ΓI) === Tuple(Γ1)
  @test size(ΓI) === (Con(D), Cov(D), Cov(D))
  @test (@inferred partials_to_Christoffel_2(∂g, I)) === ΓI

  # the first kind takes no metric at all, the second needs one
  @test_throws MethodError partials_to_Christoffel_1(∂g, g)
  @test_throws ArgumentError partials_to_Christoffel_2(∂g)
end

# Metric fields: symmetric, diagonally dominant hence invertible, with partial
# derivatives of either sign
g2F(y) = SymTensorValue{cov2}(2.0 + y[1]^2, w(1,y), 3.0 + sin(y[2]))
g3F(y) = SymTensorValue{Cov(3)}(
  2.0 + y[1]^2, w(1,y), w(2,y), 3.0 + sin(y[2]), w(3,y), 4.0 + y[3]^2)

test_Christoffel(g2F, x2)
test_Christoffel(g3F, x3)

# A metric of a dimension other than that of the partials has no index to raise
@test_throws DimensionMismatch partials_to_Christoffel_2(
  partial_derivatives(g3F,x3), inv(g2F(x2)))

# Polar coordinates, where the symbols are known in closed form: in (r,θ) with
# g = diag(1,r²) the only non-zero ones are Γʳ_θθ = -r, Γ^θ_rθ = Γ^θ_θr = 1/r,
# and Γ_rθθ = -r, Γ_θrθ = Γ_θθr = r
let r = 2.0
  g = SymTensorValue{cov2}(1.0, 0.0, r^2)
  # only ∂_r g_θθ = 2r
  ∂g = ProductTensorValue{SymNTensorValue{Tuple{cov2}},SymTensorValue{cov2}}(
    0.0, 0.0, 0.0, 0.0, 2r, 0.0)

  Γ1 = partials_to_Christoffel_1(∂g)
  Γ2 = partials_to_Christoffel_2(∂g, inv(g))

  @test Γ1[1,2,2] ≈ -r
  @test Γ1[2,1,2] ≈ r
  @test Γ1[2,2,1] ≈ r
  @test iszero(Γ1[1,1,1]) && iszero(Γ1[1,1,2]) && iszero(Γ1[2,1,1]) && iszero(Γ1[2,2,2])

  @test Γ2[1,2,2] ≈ -r
  @test Γ2[2,1,2] ≈ 1/r
  @test Γ2[2,2,1] ≈ 1/r
  @test iszero(Γ2[1,1,1]) && iszero(Γ2[1,1,2]) && iszero(Γ2[2,1,1]) && iszero(Γ2[2,2,2])
end

# A metric without independent component has constant, zero symbols
let ∂g = ProductTensorValue{SymNTensorValue{Tuple{Cov(0)}},SymTensorValue{Cov(0)},Float64}(())
  @test num_indep_components(typeof(partials_to_Christoffel_1(∂g))) == 0
  @test num_indep_components(typeof(partials_to_Christoffel_2(∂g, I))) == 0
  @test eltype(partials_to_Christoffel_1(∂g)) === Float64
end

# Integer partials widen on the halving
let ∂g = ProductTensorValue{SymNTensorValue{Tuple{cov2}},SymTensorValue{cov2},Int}(
           ntuple(n -> n, 6))
  @test eltype(partials_to_Christoffel_1(∂g)) === Float64
  @test (@inferred partials_to_Christoffel_1(∂g)) isa ProductTensorValue
end

# An order 1 first factor of either family carries the derivation index
let data = ntuple(n -> Float64(n), 6)
  ∂sym = ProductTensorValue{SymNTensorValue{Tuple{cov2}},SymTensorValue{cov2}}(data)
  ∂cov = ProductTensorValue{CoVectorValue{cov2},SymTensorValue{cov2}}(data)
  @test partials_to_Christoffel_1(∂cov) === partials_to_Christoffel_1(∂sym)
end

# The metric must be covariant, and the derivation index must share its dimension
# and variance
let data = ntuple(n -> Float64(n), 6)
  ∂up = ProductTensorValue{SymNTensorValue{Tuple{cov2}},SymTensorValue{Con(2)}}(data)
  @test_throws ArgumentError partials_to_Christoffel_1(∂up)
  @test_throws ArgumentError partials_to_Christoffel_2(∂up, inv(g2F(x2)))

  ∂con = ProductTensorValue{ContraVectorValue{Con(2)},SymTensorValue{cov2}}(data)
  @test_throws DimensionMismatch partials_to_Christoffel_1(∂con)

  ∂mix = ProductTensorValue{SymNTensorValue{Tuple{Cov(3)}},SymTensorValue{cov2}}(
           ntuple(n -> Float64(n), 9))
  @test_throws DimensionMismatch partials_to_Christoffel_1(∂mix)
end

# Codifferential in a chart, from the partials of the form, the form itself and
# the Christoffel symbols of the metric

function test_curved_codifferential(K, gF, x::Point{D}) where D
  f = y -> mk(ExteriorFormValue{K,Cov(D)}, y)
  g, ∂g = gF(x), partial_derivatives(gF, x)
  g_inv = inv(g)
  Γ = partials_to_Christoffel_2(∂g, g_inv)
  ω, ∂ω = f(x), partial_derivatives(f, x)
  δω = partials_to_codifferential(∂ω, ω, Γ, g, g_inv)

  @test typeof(δω) === ExteriorFormValue{K-1,Cov(D),Float64,binomial(D,K-1)}
  @test (@inferred partials_to_codifferential(∂ω, ω, Γ, g, g_inv)) === δω
  # the inverse metric is taken from the metric data in whichever form it comes
  @test partials_to_codifferential(∂ω, ω, Γ, g) ≈ δω
  @test partials_to_codifferential(∂ω, ω, Γ, compute_metric_data(g)) ≈ δω

  # δ is (-1)^{D(K+1)+1} ⋆d⋆ for a metric of index 0, the definition the three
  # Christoffel terms answer to: one per index of ω
  starf = y -> hodge_star(f(y), gF(y))
  dstar = partials_to_exteriorderivative(partial_derivatives(starf, x))
  @test δω ≈ (isodd(D*(K+1)+1) ? -1 : 1) * hodge_star(dstar, g)

  # vanishing symbols and the identity metric give the euclidean kernel back
  Γ0 = zero(ProductTensorValue{ContraVectorValue{Con(D)},SymTensorValue{Cov(D)},Float64})
  δI, δE = partials_to_codifferential(∂ω, ω, Γ0, I), partials_to_flat_codifferential(∂ω)
  @test typeof(δI) === typeof(δE)
  # the zero terms the full `g^{jm}` emits may carry either sign of zero
  @test δI == δE
end

for K in 1:2
  test_curved_codifferential(K, g2F, x2)
end
for K in 1:3
  test_curved_codifferential(K, g3F, x3)
end

let g = g2F(x2), ∂g = partial_derivatives(g2F,x2)
  Γ = partials_to_Christoffel_2(∂g, inv(g))

  # The codifferential of a 0-form would be a (-1)-form here too
  f0 = y -> mk(ExteriorFormValue{0,cov2}, y)
  @test_throws Exception partials_to_codifferential(
    partial_derivatives(f0,x2), f0(x2), Γ, g)

  # A form of degree above D has no component, so no metric is read at all
  fe3 = y -> zero(ExteriorFormValue{3,cov2,Float64,0})
  δe = partials_to_codifferential(partial_derivatives(fe3,x2), fe3(x2), Γ, g)
  @test typeof(δe) === ExteriorFormValue{2,cov2,Float64,1}
  @test iszero(δe)

  # The symbols have two covariant indices and one contravariant one, of the
  # dimension the form lives in
  f1 = y -> mk(CoVectorValue{cov2}, y)
  @test_throws ArgumentError partials_to_codifferential(
    partial_derivatives(f1,x2), f1(x2), partials_to_Christoffel_1(∂g), g)
  Γ3 = partials_to_Christoffel_2(partial_derivatives(g3F,x3), inv(g3F(x3)))
  @test_throws DimensionMismatch partials_to_codifferential(
    partial_derivatives(f1,x2), f1(x2), Γ3, g)
end

# Derivative rules of the metric dependent operations

# `inv`, `√|det g|` and `⋆` are not linear in the metric, so their rules read
# its value as well as its derivatives
function test_metric_partials_rules(gF, x::Point{D}) where D
  g, ∂g = gF(x), partial_derivatives(gF, x)

  inv_g, sq_g = y -> inv(gF(y)), y -> _sqrt_abs_det(gF(y))
  ∂inv = partials_op_rule(Val(1), inv)(g, ∂g)
  ∂sq  = partials_op_rule(Val(1), _sqrt_abs_det)(g, ∂g)

  @test ∂inv ≈ partial_derivatives(inv_g, x)
  @test typeof(∂inv) === typeof(partial_derivatives(inv_g, x))
  @test ∂sq ≈ partial_derivatives(sq_g, x)
  @test typeof(∂sq) === typeof(partial_derivatives(sq_g, x))

  # `⋆ω = vol ⋆ᴱ(Mᴷω)` raises the K indices of ω through the K-th compound of
  # the metric, so the rule is a product rule over `vol` and the rows of each
  # minor of that compound
  for K in 0:D, V in (Cov, Con)
    f = y -> mk(SkewSymKTensorValue{V,K,V(D)}, y)
    star = y -> hodge_star(f(y), gF(y), inv_g(y), sq_g(y))
    ∂star = partials_op_rule(Val(1), hodge_star)(
      f(x), g, inv_g(x), sq_g(x), partial_derivatives(f, x), ∂g, ∂inv, ∂sq)
    @test ∂star ≈ partial_derivatives(star, x)
    @test typeof(∂star) === typeof(partial_derivatives(star, x))
  end

  # The rule reads the metric of the variance `⋆` raises and the volume factor
  # off its arguments, and derives neither of the two from the other
  f1 = y -> mk(CoVectorValue{Cov(D)}, y)
  @test_throws ErrorException partials_op_rule(Val(1), hodge_star)(
    f1(x), g, partial_derivatives(f1, x), ∂g)
end

test_metric_partials_rules(g2F, x2)
test_metric_partials_rules(g3F, x3)

end # module
