module RotatingPΛBasesTests
# The rotating P_rΛ¹ basis (RotatingPΛBasis), built from the directional 1-form
#
#   ψ(f,k,s(α)) = dλ_k − (𝟙[k∈supp(α)]/|supp(α)|) Σ_{i∈f} dλ_i
#
# A self-contained, hand-rolled oracle (spanning set + filter + numeric
# evaluation, independent of the package internals) is implemented below and
# cross-checked against the real RotatingPΛBasis.

using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.Fields: Point, return_cache, evaluate!
using Combinatorics: combinations, multiexponents, multinomial
using LinearAlgebra
using Test

# ── Inlined oracle helpers ────────────────────────────────────────────────────

# Gridap convention: λ = (1−Σx, x…)
to_barycentric(x::Point{D,T}) where {D,T} = (one(T) - sum(x.data), x.data...)

# Project an ambient barycentric 1-form (dλ₁,…,dλ_N) to the physical D = N−1
# frame by imposing dλ₁ = −dλ₂ − … − dλ_N: phys[j] = c[j+1] − c[1].
reduce_ambient(ω::DifferentialFormValue{1,N}) where N =
    DifferentialFormValue{1,N-1}(ntuple(j -> ω.data[j+1] - ω.data[1], N-1))

# ψ(f,k,α) = dλ_k − (𝟙[k∈supp(α)] / |supp(α)|) · Σ_{i∈f} dλ_i
# returned as a coefficient vector in the AMBIENT (D+1)-dim barycentric frame
# dλ_1,…,dλ_{D+1}  (reduce_ambient later projects this to the D physical dλ's).
function rotating_psi(f::Vector{Int}, k::Int, α::NTuple{N,Int}) where N
    c = zeros(N)
    c[k] += 1.0
    supp = [i for i in 1:N if α[i] > 0]
    if !isempty(supp)
        w = (k in supp ? 1.0 : 0.0) / length(supp)
        for i in f; c[i] -= w; end
    end
    c
end

# Spanning set:  λ^α · ψ(f,k,α)
#   |α| = r EXACTLY (homogeneous barycentric degree — represents all Cartesian
#   degree-≤r polynomials via Σλᵢ=1),  f a face of Δ_D with dim(f) ≥ K,  k ∈ f,
#   supp(α) ∪ {k} = f
function rotating_spanning_set(::Val{D}, ::Val{K}, r::Int) where {D,K}
    Dp1    = D + 1
    α_list = [NTuple{Dp1,Int}(a) for a in multiexponents(Dp1, r)]
    faces  = Vector{Int}[]
    for fdim in K:D, f in combinations(1:Dp1, fdim+1)
        push!(faces, f)
    end
    spanning = Tuple{Vector{Int},Int,NTuple{Dp1,Int}}[]
    for f in faces, k in f, α in α_list
        supp = [i for i in 1:Dp1 if α[i] > 0]
        all(i -> i in f, supp) || continue
        (Set(supp) ∪ Set([k])) == Set(f) || continue
        push!(spanning, (f, k, α))
    end
    spanning
end

# Filtered basis: αᵢ = 0 for i < min(f\k)
function rotating_basis(spanning)
    filter(spanning) do (f, k, α)
        rest = [i for i in f if i != k]
        isempty(rest) && return true
        m = minimum(rest)
        all(==(0), α[i] for i in 1:(m-1))
    end
end

# Numeric evaluation of one entry at a Cartesian point x (reduced to physical D-dim 1-form)
function rotating_eval(f::Vector{Int}, k::Int, α::NTuple{N,Int}, x) where N
    λ   = to_barycentric(x)
    val = multinomial(α...) * prod(λ[i]^α[i] for i in 1:N)
    ω   = DifferentialFormValue{1,N}(Tuple(rotating_psi(f, k, α)))
    val * reduce_ambient(ω)
end

# ── Driver: D=2, K=1, r=2 ──────────────────────────────────────────────────────
D, K, r = 2, 1, 2
spanning = rotating_spanning_set(Val(D), Val(K), r)
basis    = rotating_basis(spanning)

@testset "basis size = D·C(r+D,D) = dim PrΛ¹" begin
    @test length(basis) == D * binomial(r + D, D)
end

# Rank check (numeric, reduced to physical D-dim 1-forms via reduce_ambient)
pts = [Point(0.13,0.27), Point(0.4,0.1), Point(0.05,0.6), Point(0.3,0.3),
       Point(0.6,0.2), Point(0.2,0.55), Point(0.45,0.45), Point(0.1,0.1),
       Point(0.7,0.1), Point(0.15,0.4), Point(0.5,0.05), Point(0.25,0.6)]
function rotating_rank(entries, pts)
    M = zeros(2*length(pts), length(entries))
    for (j,(f,k,α)) in enumerate(entries), (i,p) in enumerate(pts)
        v = rotating_eval(f, k, α, p)
        M[2i-1,j] = v.data[1]; M[2i,j] = v.data[2]
    end
    rank(M)
end
@testset "spanning set spans PrΛ¹ and the filtered basis is a basis" begin
    @test rotating_rank(spanning, pts) == length(basis)   # rank 𝒮 = dim PrΛ¹
    @test rotating_rank(basis, pts)    == length(basis)   # ℬ linearly independent
end

# ── Cross-check against the package implementation ────────────────────────────
#
# RotatingPΛBasis{D} is a genuine Gridap PolynomialBasis subtype generating the
# same (f,k,α) triples (via rotating_PΛ_bubbles) and the same ψ formula, but
# contracted against the simplex's actual ∂λ/∂x Jacobian instead of going
# through reduce_ambient. This test checks the two implementations agree exactly.
@testset "RotatingPΛBasis (package) == hand-rolled oracle coefficients" begin
    pkg_basis = RotatingPΛBasis(Val(D), Float64, r)
    @test length(pkg_basis) == length(basis)

    pkg_index = Dict{Tuple{Vector{Int},Int,Vector{Int}},Int}()
    for (F, bfs) in pkg_basis.bubbles, (w, k, α, _) in bfs
        pkg_index[(F, k, α)] = w
    end

    @testset "(f,k,α) triples present and ψ values match" for (f, k, α) in basis
        key = (f, k, collect(α))
        @test haskey(pkg_index, key)
        pkg_ψ    = collect(pkg_basis.Ψ[pkg_index[key]].data)
        oracle_ψ = collect(reduce_ambient(DifferentialFormValue{1,D+1}(Tuple(rotating_psi(f, k, α)))).data)
        @test pkg_ψ ≈ oracle_ψ atol=1e-12
    end
end

# ── Evaluation path: return_cache / evaluate! of the package basis ───────────
#
# The coefficient cross-check above never exercises the Gridap evaluation
# machinery (_return_cache/_setsize!/_evaluate_nd!). Evaluate the package basis
# pointwise and compare against the oracle.
@testset "RotatingPΛBasis evaluate! == oracle, pointwise" begin
    pkg_basis = RotatingPΛBasis(Val(D), Float64, r)
    pkg_index = Dict{Tuple{Vector{Int},Int,Vector{Int}},Int}()
    for (F, bfs) in pkg_basis.bubbles, (w, k, α, _) in bfs
        pkg_index[(F, k, α)] = w
    end
    cache    = return_cache(pkg_basis, pts)
    pkg_vals = evaluate!(cache, pkg_basis, pts)   # (np, ndof) matrix
    for (f, k, α) in basis, (i, p) in enumerate(pts)
        w = pkg_index[(f, k, collect(α))]
        @test collect(pkg_vals[i, w].data) ≈ collect(rotating_eval(f, k, α, p).data) atol=1e-12
    end
end

end # module
