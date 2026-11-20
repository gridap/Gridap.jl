module CovariantDerivativeTests

using Test
using Gridap.TensorValues
using LinearAlgebra: I

const cov2, cov3 = Cov(2), Cov(3)

# Deterministic components, of either sign
w(n) = sin(1.7n) + 0.3n - 0.5

mk(V) = V(ntuple(l -> w(l), num_indep_components(V)))

# Mock partial derivatives of `a`, a `ProductTensorValue` whose first factor
# carries the covariant derivation index, shaped as `partial_derivatives`
# produces. Every identity below is algebraic in `∂a` and `Γ`: they follow from
# the definition of `Γ` and from its symmetry in the two lower indices, not
# from the symmetry of mixed partials, so an arbitrary `∂a` exercises them.
∂mock(a::MultiValue, D, seed=0) =
  ProductTensorValue{SymNTensorValue{Tuple{D}},typeof(a),Float64}(
    ntuple(n -> w(n+seed), num_indep_components(typeof(a))*Int(D)))

# Metrics: symmetric and diagonally dominant, hence invertible
const g2 = SymTensorValue{cov2}(3.0, 0.4, 4.0)
const g3 = SymTensorValue{cov3}(3.0, 0.2, 0.1, 4.0, 0.3, 5.0)

christoffel(g, D, seed=0) = partials_to_Christoffel_2(∂mock(g, D, seed), inv(g))

const Γ2, Γ3 = christoffel(g2, cov2), christoffel(g3, cov3)

# ∇a_{k,O} - ∂_k a_O, one Christoffel term per index of `a`, read off the full
# tensors rather than off the components the kernel stores
function ref_covderivative_terms(a, Γ, D, k, O)
  is_cov = [S isa Cov for S in size(a)]
  sum(eachindex(O); init=zero(eltype(a))) do p
    sum(1:D; init=zero(eltype(a))) do l
      S = ntuple(i -> i == p ? l : O[i], length(O))
      (is_cov[p] ? -Γ[l,k,O[p]] : Γ[O[p],k,l]) * a[S...]
    end
  end
end

function test_covariant_derivative(a, Γ, D)
  Dc = Cov(D)
  ∂a = ∂mock(a, Dc, 3)
  u = ContraVectorValue{Con(D)}(ntuple(l -> w(l+7), D))
  ∇a  = partials_to_covariantderivative(∂a, a, Γ)
  ∇ua = partials_to_covariantderivative(u, ∂a, a, Γ)
  N = length(size(a))

  # ∇a prepends the covariant derivation index, like ∂a does, and keeps the
  # index structure of `a`, so it is a covector ⊗ K-form on a K-form
  @test size(∇a) === (Dc, size(a)...)
  @test ∇a isa ProductTensorValue
  @test size(∇ua) === size(a)
  @test (@inferred partials_to_covariantderivative(∂a, a, Γ)) === ∇a
  @test (@inferred partials_to_covariantderivative(u, ∂a, a, Γ)) === ∇ua

  for k in 1:D, O in CartesianIndices(ntuple(_->D,N))
    @test ∇a[k,Tuple(O)...] ≈
      ∂a[k,Tuple(O)...] + ref_covderivative_terms(a,Γ,D,k,Tuple(O)) atol=1e-12
  end
  # ∇_u a contracts the derivation index of ∇a with u
  for O in CartesianIndices(ntuple(_->D,N))
    @test ∇ua[Tuple(O)...] ≈ sum(u[k]*∇a[k,Tuple(O)...] for k in 1:D) atol=1e-12
  end

  # a vanishing connection leaves the partial derivatives alone
  @test partials_to_covariantderivative(∂a, a, I) === ∂a
  ∇Ia = partials_to_covariantderivative(u, ∂a, a, I)
  @test size(∇Ia) === size(a)
  for O in CartesianIndices(ntuple(_->D,N))
    @test ∇Ia[Tuple(O)...] ≈ sum(u[k]*∂a[k,Tuple(O)...] for k in 1:D) atol=1e-12
  end
end

for (D, Γ) in ((2, Γ2), (3, Γ3))
  for K in 0:D
    test_covariant_derivative(mk(ExteriorFormValue{K,Cov(D),Float64}), Γ, D)
    test_covariant_derivative(mk(ExteriorVectorValue{K,Con(D),Float64}), Γ, D)
  end
  test_covariant_derivative(mk(SymTensorValue{Cov(D),Float64}), Γ, D)
  test_covariant_derivative(mk(TensorValue{Con(D),Cov(D),Float64}), Γ, D)
end

# The Levi-Civita connection is metric compatible, ∇g = ∇g⁻¹ = 0, which pins the
# sign of each variance branch, and leaves the mixed identity alone, which pins
# them against one another. `Γ` makes these hold of the `∂g` it was built from,
# whether or not that is the derivative of anything.
function test_metric_compatibility(g, D)
  Dc = Cov(D)
  ∂g = ∂mock(g, Dc)
  Γ = partials_to_Christoffel_2(∂g, inv(g))

  @test all(c -> isapprox(c, 0.0, atol=1e-12),
            Tuple(partials_to_covariantderivative(∂g, g, Γ)))

  # the inverse metric is the other variance branch, its partials given by the
  # derivative rule of `inv` rather than by a derivative
  g_inv = inv(g)
  ∂g_inv = partials_op_rule(Val(1), inv)(g, ∂g)
  @test all(c -> isapprox(c, 0.0, atol=1e-12),
            Tuple(partials_to_covariantderivative(∂g_inv, g_inv, Γ)))

  # the mixed identity is constant, so its two Christoffel terms must cancel
  δmix = TensorValue{Con(D),Cov(D)}(ntuple(n -> Float64(fldmod1(n,D)[1] == fldmod1(n,D)[2]), D*D))
  ∂δ = zero(ProductTensorValue{SymNTensorValue{Tuple{Dc}},typeof(δmix),Float64})
  @test all(c -> isapprox(c, 0.0, atol=1e-12),
            Tuple(partials_to_covariantderivative(∂δ, δmix, Γ)))
end

test_metric_compatibility(g2, 2)
test_metric_compatibility(g3, 3)

# It is torsion free, so the Lie bracket of two vector fields is the
# antisymmetric part of their covariant derivatives: the Christoffel terms
# cancel by the symmetry of `Γ` in its lower indices
function test_torsion_free(Γ, D)
  Dc = Cov(D)
  u = ContraVectorValue{Con(D)}(ntuple(l -> w(l+1), D))
  v = ContraVectorValue{Con(D)}(ntuple(l -> w(l+13), D))
  ∂u, ∂v = ∂mock(u, Dc, 5), ∂mock(v, Dc, 17)

  ∇uv = partials_to_covariantderivative(u, ∂v, v, Γ)
  ∇vu = partials_to_covariantderivative(v, ∂u, u, Γ)
  Luv = partials_to_liederivative(u, v, ∂u, ∂v)
  @test all(isapprox.(Tuple(∇uv - ∇vu), Tuple(Luv), atol=1e-12))
end

test_torsion_free(Γ2, 2)
test_torsion_free(Γ3, 3)

# `𝑑` is the antisymmetrisation of `∇`, whose Christoffel terms cancel because
# the connection is torsion free, and `δ` its metric trace,
# `(δω)_O = -g^{jm}∇_j ω_{mO}`
function test_covariant_derivative_of_form(K, g, Γ, D)
  ω = mk(ExteriorFormValue{K,Cov(D),Float64})
  ∂ω = ∂mock(ω, Cov(D), 2)
  ∇ω = partials_to_covariantderivative(∂ω, ω, Γ)

  if K < D
    dω = partials_to_exteriorderivative(∂ω)
    for O in sorted_combinations(D, K+1)
      @test dω[O...] ≈ sum(sorting_sign(O[a], O[1:a-1]..., O[a+1:end]...) *
        ∇ω[O[a], O[1:a-1]..., O[a+1:end]...] for a in 1:K+1) atol=1e-12
    end
  end
  if K ≥ 1
    g_inv = inv(g)
    δω = partials_to_codifferential(∂ω, ω, Γ, g_inv)
    for O in sorted_combinations(D, K-1)
      @test δω[O...] ≈
        -sum(g_inv[j,m]*∇ω[j,m,O...] for j in 1:D for m in 1:D) atol=1e-12
    end
  end
end

for (D, g, Γ) in ((2, g2, Γ2), (3, g3, Γ3)), K in 0:D
  test_covariant_derivative_of_form(K, g, Γ, D)
end

# ∇ is a tensor derivation: it obeys Leibniz on a tensor product, and commutes
# with contraction. The partials of the product and of the pairing are those
# the derivative rules give, so no derivative is taken here either.
let D = 2, Γ = Γ2
  X = ContraVectorValue{Con(D)}(ntuple(l -> w(l+1), D))
  ω = mk(CoVectorValue{cov2,Float64})
  ∂X, ∂ω = ∂mock(X, cov2, 5), ∂mock(ω, cov2, 2)
  ∇X = partials_to_covariantderivative(∂X, X, Γ)
  ∇ω = partials_to_covariantderivative(∂ω, ω, Γ)

  ∂Xω = partials_op_rule(Val(1), ⊗)(X, ω, ∂X, ∂ω)
  ∇Xω = partials_to_covariantderivative(∂Xω, X⊗ω, Γ)
  for k in 1:D, i in 1:D, j in 1:D
    @test ∇Xω[k,i,j] ≈ ∇X[k,i]*ω[j] + X[i]*∇ω[k,j] atol=1e-12
  end

  # the contraction ⟨ω,X⟩ is a scalar, whose covariant derivative is ∂
  ∂pair = partials_op_rule(Val(1), ι)(X, ω, ∂X, ∂ω)
  for k in 1:D
    @test ∂pair[k] ≈ sum(∇ω[k,i]*X[i] + ω[i]*∇X[k,i] for i in 1:D) atol=1e-12
  end
end

# A scalar carries no index to transport, so ∇f is ∂f whatever the connection,
# and ∇_u f the derivative of f along u
let D = 2
  ∂f = ProductTensorValue{SymNTensorValue{Tuple{cov2}},ScalarValue{Float64},Float64}(
         ntuple(n -> w(n+4), D))
  u = ContraVectorValue{Con(D)}(0.4, -1.2)
  @test partials_to_covariantderivative(∂f, 1.5, Γ2) === ∂f
  @test partials_to_covariantderivative(∂f, 1.5, I) === ∂f
  @test partials_to_covariantderivative(u, ∂f, 1.5, Γ2) ≈ u⋅∂f
  @test partials_to_covariantderivative(u, ∂f, 1.5, I) ≈ u⋅∂f
end

# A form of degree above D has no component, and so a vanishing derivative
let D = 2, ω = zero(ExteriorFormValue{3,cov2,Float64,0})
  ∂ω = zero(ProductTensorValue{SymNTensorValue{Tuple{cov2}},typeof(ω),Float64})
  u = ContraVectorValue{Con(D)}(0.4, -1.2)
  @test num_indep_components(typeof(partials_to_covariantderivative(∂ω, ω, Γ2))) == 0
  @test typeof(partials_to_covariantderivative(u, ∂ω, ω, Γ2)) ===
    ExteriorFormValue{3,cov2,Float64,0}
end

# ∇ needs a covariant derivation index, a contravariant direction, and
# Christoffel symbols of the dimension and variance of the tensor
let ω = mk(CoVectorValue{cov2,Float64})
  ∂ω = ∂mock(ω, cov2, 2)
  u, α = ContraVectorValue{Con(2)}(0.4, -1.2), mk(CoVectorValue{cov2,Float64})

  @test_throws DimensionMismatch partials_to_covariantderivative(∂ω, ω, Γ3)
  @test_throws ArgumentError partials_to_covariantderivative(α, ∂ω, ω, Γ2)
  @test_throws ArgumentError partials_to_covariantderivative(α, ∂ω, ω, I)
  ∂con = ProductTensorValue{ContraVectorValue{Con(2)},SkewSymKTensorValue{Cov,1,cov2}}(
           ntuple(n -> Float64(n), 4))
  @test_throws ArgumentError partials_to_covariantderivative(∂con, ω, Γ2)
  @test_throws ArgumentError partials_to_covariantderivative(u, ∂con, ω, Γ2)
end

end # module
