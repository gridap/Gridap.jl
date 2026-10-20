module SerendipityPolyBasesTests

using Test
using Random
using LinearAlgebra
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials

# [Arnold & Awanou, Math. Comp. 83 (2014), (22)]
aa_dim(D,r,k) = sum(2^(D-d)*binomial(D,d)*binomial(r-d+2k,d)*binomial(d,k) for d in k:min(D,r÷2+k))
# [Gillette & Kloefkorn, Math. Comp. 88 (2019), Table 1], r = 1, 2, …
const GK_DIMS = Dict((2,1) => [4,10,17,26], (3,1) => [12,36,66,111], (3,2) => [6,21,45,82])

points(D,n) = (rng = MersenneTwister(n); [Point(ntuple(_ -> rand(rng), D)) for _ in 1:n])

# one column per function, the components at each point stacked
value_matrix(vals::AbstractMatrix) = reduce(hcat, [[c for v in col for c in Tuple(v)] for col in eachcol(vals)])
value_matrix(b::AbstractVector{<:Field}, x) = value_matrix(evaluate(b, x))
value_matrix(fs::AbstractVector{<:Function}, x) = value_matrix([f(xi) for xi in x, f in fs])

in_span(A,B) = rank([A B]; rtol=1e-10) == rank(A; rtol=1e-10)
same_span(A,B) = in_span(A,B) && in_span(B,A)

monomials(D,r) = [x -> prod(Tuple(x).^α) for α in Iterators.product(ntuple(_ -> 0:r, D)...) if sum(α) <= r]
vector_monomials(D,r) = [x -> m(x)*e for m in monomials(D,r) for e in component_basis(VectorValue{D,Float64})]
rot(v) = VectorValue(-v[2], v[1])

T = Float64

@test Polynomials._default_poly_type(:S⁻,0,2) == ModalC0
@test Polynomials._default_poly_type(:S,1,2) == Legendre
@test Polynomials._default_poly_type(:S⁻,2,3) == Legendre
@test Polynomials._default_poly_type(:S,0,2) == ModalC0
@test Polynomials._default_poly_type(:P,1,2) == Bernstein

# 𝓢⁻ᵣΛ⁰ = 𝓢ᵣΛ⁰ and 𝓢⁻ᵣΛᴰ = 𝓟ᵣ₋₁Λᴰ
for D in 1:3
  x = points(D,10)
  @test evaluate(FEEC_poly_basis(Val(D),T,3,0,:S⁻,Monomial),x) == evaluate(FEEC_poly_basis(Val(D),T,3,0,:S,Monomial),x)
  @test evaluate(FEEC_poly_basis(Val(D),T,3,D,:S⁻,Monomial),x) == evaluate(MonomialBasis(Val(D),T,2,Polynomials._p_filter),x)
end

# Lengths, types, independence
for (D,k) in ((2,1),(3,1),(3,2)), r in 1:(D == 2 ? 4 : 3)
  V = VectorValue{D,T}
  x = points(D, 3*aa_dim(D,r,k))

  b = FEEC_poly_basis(Val(D),T,r,k,:S)
  @test b isa SerendipityPolyBasis{D,V,Legendre}
  @test length(b) == aa_dim(D,r,k)
  @test rank(value_matrix(b,x)) == length(b)

  bt = FEEC_poly_basis(Val(D),T,r,k,:S⁻)
  @test bt isa SerendipityPolyBasis{D,V,Legendre}
  @test length(bt) == GK_DIMS[(D,k)][r]
  @test rank(value_matrix(bt,x)) == length(bt)

  # the index sets define the space for any hierarchical family
  for PT in (Monomial, Chebyshev{:T})
    @test same_span(value_matrix(b,x), value_matrix(FEEC_poly_basis(Val(D),T,r,k,:S,PT),x))
    @test same_span(value_matrix(bt,x), value_matrix(FEEC_poly_basis(Val(D),T,r,k,:S⁻,PT),x))
  end
end

b = SerendipityPolyBasis{3}(Legendre,T,2,1)
@test testvalue(typeof(b)) isa typeof(b)
@test get_order(b) == 4
@test get_orders(b) == (4,4,4)
x = points(3,5)
@test_throws ErrorException evaluate(Broadcasting(∇∇)(b),x)
@test_throws AssertionError SerendipityPolyBasis{3}(Bernstein,T,2,1)
@test_throws AssertionError SerendipityPolyBasis{3}(Legendre,T,2,3)

# Gradients
for (D,k,trimmed) in ((2,1,false),(2,1,true),(3,1,false),(3,2,true))
  b = SerendipityPolyBasis{D}(Legendre,T,3,k; trimmed)
  x = points(D,4)
  ∇b = evaluate(Broadcasting(∇)(b),x)
  h = 1e-6
  for q in 1:D
    e = Point(ntuple(d -> d == q ? h : 0.0, D))
    fd = (evaluate(b,x .+ Ref(e)) - evaluate(b,x .- Ref(e)))/(2h)
    @test maximum(abs, map((g,v) -> maximum(abs, Tuple(g)[q:D:end] .- Tuple(v)), ∇b, fd)) < 1e-6
  end
end

# The spaces [Arnold & Awanou, §1]: in 2D, 𝓟ᵣ(□;R²) ⊕ d{x^{r+1}y, xy^{r+1}}, with
# either proxy; in 3D, by hand, 𝓟ᵣ(□;R³) plus
#   H(div):  curl(x₂x₃(w₂-w₃), x₃x₁(w₃-w₁), x₁x₂(w₁-w₂)), wᵢ ∈ 𝓟ᵣ, independent of xᵢ
#   H(curl): (x₂x₃(w₂-w₃), x₃x₁(w₃-w₁), x₁x₂(w₁-w₂)), wᵢ ∈ 𝓟ᵣ₋₁, independent of xᵢ,
#            and ∇s for s of superlinear degree ≤ r+1
function aa_potentials(r)
  map(Iterators.product(1:3, monomials(2,r))) do (i,m)
    w(x) = m(Point(Tuple(x)[[j for j in 1:3 if j != i]]...))
    x -> VectorValue(x[2]*x[3]*((i == 2) - (i == 3))*w(x),
                     x[3]*x[1]*((i == 3) - (i == 1))*w(x),
                     x[1]*x[2]*((i == 1) - (i == 2))*w(x))
  end |> vec
end
sdeg(α) = sum(a for a in α if a >= 2; init=0)

for r in 1:4
  x = points(2,60)
  grads = [∇(x -> x[1]^(r+1)*x[2]), ∇(x -> x[1]*x[2]^(r+1))]
  aa_curl = [vector_monomials(2,r); grads]
  aa_div = [vector_monomials(2,r); [x -> rot(g(x)) for g in grads]]
  @test same_span(value_matrix(FEEC_poly_basis(Val(2),T,r,1,:S),x), value_matrix(aa_curl,x))
  @test same_span(value_matrix(FEEC_poly_basis(Val(2),T,r,1,:S; rotate_90=true),x), value_matrix(aa_div,x))
end
for r in 1:3
  x = points(3,3*aa_dim(3,r,1))
  aa_div = [vector_monomials(3,r); [∇×W for W in aa_potentials(r)]]
  @test same_span(value_matrix(FEEC_poly_basis(Val(3),T,r,2,:S),x), value_matrix(aa_div,x))
  grads = [∇(x -> prod(Tuple(x).^α)) for α in Iterators.product(ntuple(_ -> 0:(r+3), 3)...)
           if sdeg(α) <= r+1 && sum(α) > r+1]
  aa_curl = [vector_monomials(3,r); aa_potentials(r-1); grads]
  @test same_span(value_matrix(FEEC_poly_basis(Val(3),T,r,1,:S),x), value_matrix(aa_curl,x))
end

# The trimmed spaces [Gillette & Kloefkorn, (3.1)]: 𝓢⁻ᵣΛᵏ = 𝓢ᵣ₋₁Λᵏ + κ𝓢ᵣ₋₁Λᵏ⁺¹,
# κ being x⋅ on 𝓢ᵣ₋₁Λᴰ = 𝓟ᵣ₋₁ and x× on 𝓢ᵣ₋₁Λ² in 3D; 𝓠⁻₁Λᵏ at r = 1
for (D,k) in ((2,1),(3,1),(3,2)), r in 2:3
  x = points(D,3*aa_dim(D,r,k))
  S = evaluate(FEEC_poly_basis(Val(D),T,r-1,k,:S; rotate_90=(D == 2)),x)
  if k == D-1
    κS = value_matrix([y -> m(y)*y for m in monomials(D,r-1)],x)
  else
    w = evaluate(FEEC_poly_basis(Val(D),T,r-1,k+1,:S),x)
    κS = value_matrix([cross(x[q],w[q,j]) for q in axes(w,1), j in axes(w,2)])
  end
  bt = FEEC_poly_basis(Val(D),T,r,k,:S⁻; rotate_90=(D == 2))
  @test same_span(value_matrix(bt,x), [value_matrix(S) κS])
end
for (D,k) in ((2,1),(3,1),(3,2))
  x = points(D,30)
  @test same_span(value_matrix(FEEC_poly_basis(Val(D),T,1,k,:S⁻),x),
                  value_matrix(FEEC_poly_basis(Val(D),T,1,k,:Q⁻),x))
end

# The complex: ∇𝓢ᵣ₊₁Λ⁰ ⊂ 𝓢ᵣΛ¹, curl 𝓢ᵣΛ¹ ⊂ 𝓢ᵣ₋₁Λᴰ⁻¹, div 𝓢ᵣΛᴰ⁻¹ = 𝓟ᵣ₋₁,
# and likewise ∇𝓢ᵣΛ⁰ ⊂ 𝓢⁻ᵣΛ¹, curl 𝓢⁻ᵣΛ¹ ⊂ 𝓢⁻ᵣΛᴰ⁻¹, div 𝓢⁻ᵣΛᴰ⁻¹ = 𝓟ᵣ₋₁
for D in (2,3), r in 1:3, F in (:S,:S⁻)
  x = points(D,150)
  r0 = F == :S ? r+1 : r
  curl_b = FEEC_poly_basis(Val(D),T,r,1,F)
  div_b = FEEC_poly_basis(Val(D),T,r,D-1,F; rotate_90=(D == 2))
  ∇S0 = evaluate(Broadcasting(∇)(FEEC_poly_basis(Val(D),T,r0,0,:S)),x)
  @test in_span(value_matrix(curl_b,x), value_matrix(∇S0))
  @test same_span(value_matrix(evaluate(Broadcasting(divergence)(div_b),x)), value_matrix(monomials(D,r-1),x))
  curls = value_matrix(evaluate(Broadcasting(curl)(curl_b),x))
  if D == 2
    @test same_span(curls, value_matrix(monomials(2,r-1),x))
  elseif F == :S⁻ || r > 1
    target = FEEC_poly_basis(Val(D),T,ifelse(F == :S, r-1, r),D-1,F)
    @test in_span(value_matrix(target,x), curls)
  end
end

# Exact forms d𝓗ₛ₊₁Λᵏ⁻¹, modulo 𝓟ₛ₋₁Λᵏ
dimH(D,s) = binomial(s+D-1,D-1)
for s in 0:3
  x = points(3,60)
  g = Polynomials._exact_poly_basis(Legendre,T,3,s,1)     # ∇𝓗ₛ₊₁
  @test length(g) == dimH(3,s+1)
  @test rank(value_matrix(g,x)) == length(g)
  @test maximum(abs, value_matrix(evaluate(Broadcasting(curl)(g),x))) < 1e-10
  c = Polynomials._exact_poly_basis(Legendre,T,3,s,2)     # curl 𝓗ₛ₊₁(R³)
  @test length(c) == 3*dimH(3,s+1) - dimH(3,s+2) + (s == -1)
  @test rank(value_matrix(c,x)) == length(c)
  @test maximum(abs, evaluate(Broadcasting(divergence)(c),x)) < 1e-10
  r = Polynomials._exact_poly_basis(Legendre,T,2,s,1; rotate_90=true) # curl 𝓗ₛ₊₁(R²)
  @test length(r) == dimH(2,s+1)
  @test maximum(abs, evaluate(Broadcasting(divergence)(r),points(2,20))) < 1e-10
end

end # module
