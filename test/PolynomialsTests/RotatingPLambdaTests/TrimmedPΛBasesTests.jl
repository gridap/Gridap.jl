module TrimmedPΛBasesTests
# The trimmed P_r⁻Λ¹ basis, built from the Whitney (directional) 1-form
#
#   ϕ(ξ;e1,e2) = ξ^{e1} dξ^{e2} − ξ^{e2} dξ^{e1}
#
# A self-contained, hand-rolled oracle (spanning set + filter + numeric
# evaluation, independent of the package internals) is implemented below and
# cross-checked against BarycentricPmΛBasis with flavor=:BMM. Unlike ψ in the
# untrimmed case (RotatingPΛBasesTests.jl), ϕ is NOT constant-coefficient — it
# is linear in λ — so the cross-check compares pointwise numeric evaluations
# rather than constant coefficient vectors.

using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.Fields: Point, return_cache, evaluate!
using Combinatorics: combinations, multiexponents, multinomial
using LinearAlgebra
using Test

# ── Inlined oracle helpers ────────────────────────────────────────────────────

# Gridap convention: λ = (1−Σx, x…)
to_barycentric(x::Point{D,T}) where {D,T} = (one(T) - sum(x.data), x.data...)

# Project an ambient barycentric 1-form (dλ¹,…,dλᴺ) to the physical D = N−1
# frame by imposing dλ¹ = −dλ² − … − dλᴺ: phys[j] = c[j+1] − c[1].
reduce_ambient(ω::ExteriorFormValue{1,N}) where N =
    ExteriorFormValue{1,N-1}(ntuple(j -> ω.data[j+1] - ω.data[1], N-1))

# ϕ(e1,e2;λ) = λ^{e1} dλ^{e2} − λ^{e2} dλ^{e1}, returned as a coefficient vector in
# the AMBIENT (D+1)-dim barycentric frame dλ¹,…,dλ^{D+1} at the barycentric
# point λ (length N). Position-dependent, unlike rotating_psi.
function trimmed_phi(e1::Int, e2::Int, N::Int, λ)
    c = zeros(eltype(λ), N)
    c[e2] += λ[e1]
    c[e1] -= λ[e2]
    c
end

# Spanning set:  λ^α · ϕ(e1,e2)
#   f a face of Δ_D with dim(f) ≥ K (K=1: a pair needs ≥2 vertices), e=(e1,e2)
#   a pair of DISTINCT vertices of f (both orders included), |α| = r-1 EXACTLY,
#   supp(α) ∪ {e1,e2} = f
function trimmed_spanning_set(::Val{D}, ::Val{K}, r::Int) where {D,K}
    Dp1    = D + 1
    α_list = [NTuple{Dp1,Int}(a) for a in multiexponents(Dp1, r-1)]
    faces  = Vector{Int}[]
    for fdim in K:D, f in combinations(1:Dp1, fdim+1)
        push!(faces, f)
    end
    spanning = Tuple{Vector{Int},Tuple{Int,Int},NTuple{Dp1,Int}}[]
    for f in faces, e1 in f, e2 in f, α in α_list
        e1 == e2 && continue
        supp = [i for i in 1:Dp1 if α[i] > 0]
        all(i -> i in f, supp) || continue
        (Set(supp) ∪ Set([e1,e2])) == Set(f) || continue
        push!(spanning, (f, (e1,e2), α))
    end
    spanning
end

# Filtered basis: e1 < e2 and e1 = min(f)
function trimmed_basis(spanning)
    filter(spanning) do (f, e, α)
        e1, e2 = e
        e1 < e2 && e1 == minimum(f)
    end
end

# Numeric evaluation of one entry at a Cartesian point x (reduced to physical D-dim 1-form)
function trimmed_eval(f::Vector{Int}, e::Tuple{Int,Int}, α::NTuple{N,Int}, x) where N
    λ      = to_barycentric(x)
    e1, e2 = e
    val    = prod(λ[i]^α[i] for i in 1:N)   # bare monomial, i.e. flavor=:BMM
    ω      = ExteriorFormValue{1,N}(Tuple(trimmed_phi(e1, e2, N, λ)))
    val * reduce_ambient(ω)
end

# ── Driver: D=2, K=1, r=3 ──────────────────────────────────────────────────────
# r must be ≥ 3 for the bare monomial λ^α of flavor=:BMM to differ from the
# Bernstein B_α of :AFW: below that |α| = r-1 ≤ 1 and multinomial(α) is 1.
D, K, r = 2, 1, 3
spanning = trimmed_spanning_set(Val(D), Val(K), r)
basis    = trimmed_basis(spanning)

@testset "basis size = r·C(r+D,D−1) = dim Pr⁻Λ¹" begin
    @test length(basis) == r * binomial(r + D, D - 1)
end

# Rank check (numeric, reduced to physical D-dim 1-forms via reduce_ambient)
pts = [Point(0.13,0.27), Point(0.4,0.1), Point(0.05,0.6), Point(0.3,0.3),
       Point(0.6,0.2), Point(0.2,0.55), Point(0.45,0.45), Point(0.1,0.1),
       Point(0.7,0.1), Point(0.15,0.4), Point(0.5,0.05), Point(0.25,0.6)]
function trimmed_rank(entries, pts)
    M = zeros(2*length(pts), length(entries))
    for (j,(f,e,α)) in enumerate(entries), (i,p) in enumerate(pts)
        v = trimmed_eval(f, e, α, p)
        M[2i-1,j] = v.data[1]; M[2i,j] = v.data[2]
    end
    rank(M)
end
@testset "spanning set spans Pr⁻Λ¹ and the filtered basis is a basis" begin
    @test trimmed_rank(spanning, pts) == length(basis)   # rank 𝒮⁻ = dim Pr⁻Λ¹
    @test trimmed_rank(basis, pts)    == length(basis)   # ℬ⁻ linearly independent
end

# ── Cross-check against the package implementation ────────────────────────────
#
# BarycentricPmΛBasis{D} generates the same (f,e,α) triples — its bubble
# functions carry the pair as J, and its filter α_i = 0 for i < min(J) forces
# min(J) = min(f), which is the oracle's e1 = min(f) — and the same ϕ formula,
# but contracted against the simplex's actual ∂λ/∂x Jacobian instead of going
# through reduce_ambient. Since ϕ is position-dependent (unlike ψ), the two
# implementations are compared pointwise at several sample Cartesian points.
#
# flavor=:BMM is what scales by the bare monomial λ^α of the oracle; :AFW
# scales by B_α instead, and the last testset pins that difference down.
pkg_index(b) = Dict((F, (J[1], J[2]), α) => w
                    for (F, bfs) in get_bubbles(b) for (w, α, _, J) in bfs)

@testset "BarycentricPmΛBasis(:BMM) == hand-rolled oracle, pointwise" begin
    pkg_basis = BarycentricPmΛBasis(Val(D), Float64, r, K; flavor=:BMM)
    @test length(pkg_basis) == length(basis)

    index    = pkg_index(pkg_basis)
    cache    = return_cache(pkg_basis, pts)
    pkg_vals = evaluate!(cache, pkg_basis, pts)   # (np, ndof) matrix of VectorValue{D}

    @testset "(f,e,α) triples present and ϕ values match pointwise" for (f, e, α) in basis
        key = (f, e, collect(α))
        @test haskey(index, key)
        w = index[key]
        for (i, p) in enumerate(pts)
            pkg_v    = collect(pkg_vals[i,w].data)
            oracle_v = collect(trimmed_eval(f, e, α, p).data)
            @test pkg_v ≈ oracle_v atol=1e-12
        end
    end
end

@testset "flavor=:AFW scales the same ϕ by B_α instead of λ^α" begin
    afw = BarycentricPmΛBasis(Val(D), Float64, r, K; flavor=:AFW)
    index = pkg_index(afw)
    afw_vals = evaluate!(return_cache(afw, pts), afw, pts)

    hits = 0
    for (f, e, α) in basis
        w = index[(f, e, collect(α))]
        c = multinomial(α...)
        c == 1 && continue
        hits += 1
        for (i, p) in enumerate(pts)
            @test collect(afw_vals[i,w].data) ≈ c .* collect(trimmed_eval(f, e, α, p).data) atol=1e-12
        end
    end
    @test hits > 0   # otherwise the two flavors coincide and this proves nothing
end

end # module
