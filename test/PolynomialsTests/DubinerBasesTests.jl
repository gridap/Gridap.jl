module DubinerBasesTests

using LinearAlgebra
using Random: MersenneTwister
using Test

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs

############################################################################################
# Helpers
############################################################################################

# The Dubiner polynomial of multi-index α, written the textbook way: in
# collapsed ("Duffy") coordinates ηₖ = xₖ/σ_{k-1}, with an explicit Jacobi
# polynomial per factor. This is the definition `DubinerBases.jl` implements, and
# it is *not* how it implements it -- the σ^α factors and the division by σ are
# replaced there by a homogeneous recurrence. Valid only where every σ > 0,
# which is why the vertex tests below use the monomial fit instead.
function dubiner_reference(α::NTuple{D,Int}, x) where D
  v, A = 1.0, 0
  for k in 1:D
    σ = 1 - sum(x[i] for i in (k + 1):D; init=0.0)   # σ_k = x_k + σ_{k-1}
    a = 2*A + (k - 1)
    η = x[k] / σ
    v *= σ^α[k] * Polynomials.jacobi(2η - 1, α[k], Float64(a), 0.0)
    A += α[k]
  end
  v * Polynomials._dubiner_scale(α)
end

# Points strictly inside the reference D-simplex.
function interior_points(::Val{D}, n; seed=1) where D
  rng = MersenneTwister(seed)
  pts = Point{D,Float64}[]
  while length(pts) < n
    y = ntuple(_ -> rand(rng), Val(D))
    sum(y) < 0.98 && push!(pts, Point(y))
  end
  pts
end

# The matrix expressing a basis in the monomials of P_K, by an (exact,
# overdetermined) fit; and the fit's residual, which must be at round-off if the
# basis really is in P_K.
function monomial_fit(b, ::Val{D}, K, pts) where D
  mb = MonomialBasis(Val(D), Float64, K, Polynomials._p_filter)
  M, B = evaluate(mb, pts), evaluate(b, pts)
  C = M \ B
  (mb, C, maximum(abs, M * C - B))
end

function gram(b, p, degree)
  q = Quadrature(p, degree)
  V = evaluate(b, get_coordinates(q))
  transpose(V) * (get_weights(q) .* V)
end

############################################################################################
# Basis API
############################################################################################

for (D, K) in ((1, 3), (2, 4), (3, 3))
  b = DubinerBasis(Val(D), Float64, K)
  @test length(b) == binomial(D + K, D)
  @test get_order(b) == K
  @test get_orders(b) == ntuple(_ -> K, D)
  @test value_type(b) == Float64
  @test get_dimension(b) == D
  @test all(α -> sum(α) <= K, get_exponents(b))
  @test length(unique(get_exponents(b))) == length(b)
end

@test isHierarchical(Dubiner)
@test length(testvalue(DubinerBasis{2,Float64})) == 1

# a MultiValue basis is the direct sum of one scalar copy per component
for V in (VectorValue{2,Float64}, TensorValue{2,2,Float64,4},
          SymTensorValue{2,Float64,3})
  b = DubinerBasis(Val(2), V, 2)
  @test length(b) == 6*num_indep_components(V)
  @test value_type(b) == V
end

############################################################################################
# The whole construction collapses to Gridap's own basis when D = 1, which pins
# both the normalisation and the [0,1] reference interval
############################################################################################

for K in 0:5
  b = DubinerBasis(Val(1), Float64, K)
  l = LegendreBasis(Val(1), Float64, K)
  x = [Point(t) for t in range(0, 1; length=11)]
  @test evaluate(b, x) ≈ evaluate(l, x)
  @test evaluate(Broadcasting(∇)(b), x) ≈ evaluate(Broadcasting(∇)(l), x)
end

############################################################################################
# Orthonormality
############################################################################################

for K in 0:5
  @test norm(gram(DubinerBasis(Val(2), Float64, K), TRI, 2*K + 2) - I) < 1e-13
end
for K in 0:4
  @test norm(gram(DubinerBasis(Val(3), Float64, K), TET, 2*K + 2) - I) < 1e-13
end

# the mean of every non-constant member vanishes: the property the whole
# exercise is about, and what `LegendreBasis` does NOT give on a simplex
for (p, D) in ((TRI, 2), (TET, 3))
  b = DubinerBasis(Val(D), Float64, 3)
  q = Quadrature(p, 8)
  means = vec(transpose(get_weights(q)) * evaluate(b, get_coordinates(q)))
  @test all(i -> abs(means[i]) < 1e-14, 2:length(b))
  lb = LegendreBasis(Val(D), Float64, 3, Polynomials._p_filter)
  lmeans = vec(transpose(get_weights(q)) * evaluate(lb, get_coordinates(q)))
  @test maximum(abs, lmeans[2:end]) > 0.1     # ... unlike Legendre
end

############################################################################################
# Against the collapsed-coordinate definition
#
# An independent implementation, straight from the textbook formula with
# `PolynomialBases.jacobi`, checking the homogeneous recurrence.
############################################################################################

for (D, K) in ((2, 5), (3, 4))
  b = DubinerBasis(Val(D), Float64, K)
  pts = interior_points(Val(D), 20)
  ref = [dubiner_reference(α, x) for x in pts, α in get_exponents(b)]
  @test evaluate(b, pts) ≈ ref atol = 1e-11
end

############################################################################################
# At the collapsed vertex
#
# ηₖ = xₖ/σ_{k-1} is 0/0 at the vertex eD; the polynomials are perfectly regular
# there and the homogeneous form must return their value.
############################################################################################

for (p, D, K) in ((TRI, 2, 5), (TET, 3, 4))
  b = DubinerBasis(Val(D), Float64, K)
  verts = get_vertex_coordinates(p)
  @test all(isfinite, evaluate(b, verts))

  # against the monomial expansion, which knows nothing of collapsed coordinates
  pts = interior_points(Val(D), 200)
  _, C, res = monomial_fit(b, Val(D), K, pts)
  @test res < 1e-10
  mb = MonomialBasis(Val(D), Float64, K, Polynomials._p_filter)
  @test evaluate(b, verts) ≈ evaluate(mb, verts) * C atol = 1e-9
end

############################################################################################
# The basis spans P_K
############################################################################################

for (D, K) in ((2, 5), (3, 4))
  b = DubinerBasis(Val(D), Float64, K)
  pts = interior_points(Val(D), 300)
  _, C, res = monomial_fit(b, Val(D), K, pts)
  @test res < 1e-10
  @test size(C) == (binomial(D + K, D), length(b))
  @test rank(C; atol=1e-8) == length(b)
end

############################################################################################
# Gradients
#
# Differentiate the monomial expansion instead of the basis: an independent route
# to ∇φ that shares no code with the differentiated recurrence.
############################################################################################

for (D, K) in ((2, 5), (3, 4))
  b = DubinerBasis(Val(D), Float64, K)
  pts = interior_points(Val(D), 200)
  mb, C, res = monomial_fit(b, Val(D), K, pts)
  @test res < 1e-10
  gb = evaluate(Broadcasting(∇)(b), pts)
  gm = evaluate(Broadcasting(∇)(mb), pts)
  err = maximum(norm(gb[i, j] - sum(C[m, j] * gm[i, m] for m in axes(C, 1)))
                for i in axes(gb, 1), j in axes(gb, 2))
  @test err < 1e-9
end

# gradients at the collapsed vertex too
b4 = DubinerBasis(Val(2), Float64, 4)
@test all(isfinite, norm.(evaluate(Broadcasting(∇)(b4), get_vertex_coordinates(TRI))))
@test_throws ErrorException evaluate(Broadcasting(∇∇)(b4), [Point(0.2, 0.3)])

############################################################################################
# MultiValue
#
# `_cartprod_set_value!` / `_cartprod_set_derivative!` are reused verbatim, so
# this only checks that a MultiValue basis is the scalar one scattered.
############################################################################################

mv_pts = interior_points(Val(2), 12)
sb = DubinerBasis(Val(2), Float64, 3)
sv = evaluate(sb, mv_pts)
for V in (VectorValue{2,Float64}, SymTensorValue{2,Float64,3})
  n = num_indep_components(V)
  vb = DubinerBasis(Val(2), V, 3)
  vv = evaluate(vb, mv_pts)
  @test size(vv) == (length(mv_pts), length(sb) * n)
  for j in axes(sv, 2), c in 1:n
    comp = [vv[i, n * (j - 1) + c] for i in axes(sv, 1)]
    expect = [V(ntuple(q -> q == c ? sv[i, j] : 0.0, n)...) for i in axes(sv, 1)]
    @test comp ≈ expect
  end
end

############################################################################################
# Filters
#
# The payoff: with an orthogonal basis, `Pᵏ ∩ (Pᵠ)^⊥` is a *selection*.
############################################################################################

q2 = Quadrature(TRI, 10)
xq2, wq2 = get_coordinates(q2), get_weights(q2)
for (k, qq) in ((2, 0), (2, 1), (3, 1), (4, 2))
  b = DubinerBasis(Val(2), Float64, k, Polynomials._pk_minus_pq_filter(k, qq))
  @test length(b) == binomial(2 + k, 2) - binomial(2 + qq, 2)
  lower = MonomialBasis(Val(2), Float64, qq, Polynomials._p_filter)
  M = transpose(evaluate(lower, xq2)) * (wq2 .* evaluate(b, xq2))
  @test maximum(abs, M) < 1e-14        # orthogonal to P_q, term by term
  @test norm(gram(b, TRI, 2*k + 2) - I) < 1e-13
end

# in 3D as well
b3 = DubinerBasis(Val(3), Float64, 2, Polynomials._pk_minus_pq_filter(2, 0))
q3 = Quadrature(TET, 8)
means3 = vec(transpose(get_weights(q3)) * evaluate(b3, get_coordinates(q3)))
@test length(b3) == 10 - 1
@test maximum(abs, means3) < 1e-14

############################################################################################
# The filtered basis reproduces hand-built cell weights
#
# A basis of P₂ ∩ P₀^⊥ can be built by subtracting each monomial's mean, and one
# of P₂ ∩ P₁^⊥ by L² projection. Both constructions are reproduced here -- the
# span is all that enters a constraint functional, so this is the statement that
# selecting from the Dubiner basis gives the same space.
############################################################################################

qw = Quadrature(TRI, 10)
xq, wq = get_coordinates(qw), get_weights(qw)
mono = evaluate(MonomialBasis(Val(2), Float64, 2, Polynomials._p_filter), xq)

# P₂ ∩ P₀^⊥, by mean subtraction
mono_means = vec(transpose(wq) * mono) ./ sum(wq)
hand0 = hcat((mono[:, i] .- mono_means[i] .* mono[:, 1] for i in 2:6)...)

# P₂ ∩ P₁^⊥, by L² projection onto the complement of P₁
lin = evaluate(MonomialBasis(Val(2), Float64, 1, Polynomials._p_filter), xq)
R = mono - lin * ((transpose(lin) * (wq .* lin)) \ (transpose(lin) * (wq .* mono)))
F = eigen(Symmetric(transpose(R) * (wq .* R)))
hand1 = R * F.vectors[:, findall(λ -> λ > 1e-10*maximum(F.values), F.values)]

for (qq, A) in ((0, hand0), (1, hand1))
  B = evaluate(DubinerBasis(Val(2), Float64, 2, Polynomials._pk_minus_pq_filter(2, qq)), xq)
  @test size(A, 2) == size(B, 2)
  # same span: each is an exact linear combination of the other
  @test maximum(abs, A * (A \ B) - B) < 1e-11
  @test maximum(abs, B * (B \ A) - A) < 1e-11
  # and both are orthogonal to P_q
  L = evaluate(MonomialBasis(Val(2), Float64, qq, Polynomials._p_filter), xq)
  @test maximum(abs, transpose(L) * (wq .* B)) < 1e-14
  @test maximum(abs, transpose(L) * (wq .* A)) < 1e-13
end

############################################################################################
# Conditioning
#
# The reason to want an orthogonal basis beyond tidiness: the generalized
# Vandermonde against point evaluations stays far better conditioned than the
# monomial one as the degree grows.
############################################################################################

for K in (3, 5, 7)
  b = DubinerBasis(Val(2), Float64, K)
  m = MonomialBasis(Val(2), Float64, K, Polynomials._p_filter)
  nodes = interior_points(Val(2), binomial(2 + K, 2); seed=7)
  @test cond(evaluate(b, nodes)) < cond(evaluate(m, nodes))
end

end # module
