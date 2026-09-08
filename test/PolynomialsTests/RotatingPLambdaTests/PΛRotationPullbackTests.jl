module PΛRotationPullbackTests
# End-to-end validation of the rotation change of basis against a NUMERIC
# pullback, for the untrimmed BarycentricPΛBasis 1-forms and the trimmed
# BarycentricPmΛBasis ones, in both flavors each.
#
# The self-consistency tests in PΛRotationsTests.jl / PΛTrimmedRotationsTests.jl
# check the closed-form index calculus against itself; this file closes the
# loop with actual function evaluation. For the affine simplex automorphism
# A_p (mapping vertex V_i ↦ V_{p[i]}) with constant Jacobian J_p, the pullback
# of a 1-form field w is (A_p^* w)(x) = J_p' · w(A_p x), and the covariance of
# the spanning family is A_p^* w(ξ;t) = w(ξ; p⁻¹·t). Since the closed form
# targets the π-pushforward entries, C = rotation_change_of_basis(b, π) is the
# matrix of the pullback along A_{π⁻¹} expressed in the FIXED basis:
#
#   (A_{invperm(π)}^* w_μ)(x) = Σ_ν C[μ,ν] · w_ν(x)      (no relabelling)
#
# The test evaluates both sides numerically (self-contained oracle formulas,
# independent of the basis implementation), asserts the identity for ALL π
# (both spaces, several (D,r)), and asserts the DIRECTION: with the map A_π
# instead of A_{π⁻¹} the identity must fail for every π with π² ≠ id (a
# wrong-but-self-consistent C would fail both).

using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.Polynomials: bubble_entries, rotate_basis_function, rotation_change_of_basis
using Gridap.Polynomials: RotationCache, rotation_map
using Gridap.Fields: Point
using Combinatorics: permutations, multinomial
using LinearAlgebra
using Test

# ── Inlined helpers (from the ExteriorGridap.jl SimplexFEECBases oracle) ──────

# Gridap convention: λ = (1−Σx, x…)
to_barycentric(x::Point{D,T}) where {D,T} = (one(T) - sum(x.data), x.data...)

# Project an ambient barycentric 1-form (dλ¹,…,dλᴺ) to the physical D = N−1
# frame by imposing dλ¹ = −dλ² − … − dλᴺ: phys[j] = c[j+1] − c[1].
reduce_ambient(ω::ExteriorFormValue{1,N}) where N =
    ExteriorFormValue{1,N-1}(ntuple(j -> ω.data[j+1] - ω.data[1], N-1))

# ── Oracles (identical formulas to the basis test files) ──────────────────────

function rotating_psi(f::Vector{Int}, k::Int, α::Vector{Int}, N::Int)
    c = zeros(N)
    c[k] += 1.0
    supp = [i for i in 1:N if α[i] > 0]
    if !isempty(supp)
        w = (k in supp ? 1.0 : 0.0) / length(supp)
        for i in f; c[i] -= w; end
    end
    c
end

function oracle_eval(f::Vector{Int}, k::Int, α::Vector{Int}, x)   # untrimmed
    λ = to_barycentric(x); N = length(λ)
    val = multinomial(α...) * prod(λ[i]^α[i] for i in 1:N)
    ω = ExteriorFormValue{1,N}(Tuple(rotating_psi(f, k, α, N)))
    val * reduce_ambient(ω)
end

# :AFW direction form: α itself weights the correction, instead of its support.
function afw_phi(f::Vector{Int}, k::Int, α::Vector{Int}, N::Int)
    c = zeros(N)
    c[k] += 1.0
    r = sum(α)
    for i in f; c[i] -= α[k]/r; end
    c
end

function oracle_eval_afw(f::Vector{Int}, k::Int, α::Vector{Int}, x)
    λ = to_barycentric(x); N = length(λ)
    val = multinomial(α...) * prod(λ[i]^α[i] for i in 1:N)
    ω = ExteriorFormValue{1,N}(Tuple(afw_phi(f, k, α, N)))
    val * reduce_ambient(ω)
end

function oracle_eval(f::Vector{Int}, e::Tuple{Int,Int}, α::Vector{Int}, x)   # trimmed
    λ = to_barycentric(x); N = length(λ)
    e1, e2 = e
    c = zeros(eltype(λ), N); c[e2] += λ[e1]; c[e1] -= λ[e2]
    val = prod(λ[i]^α[i] for i in 1:N)   # BARE monomial λ^α, the :BMM scaling
    val * reduce_ambient(ExteriorFormValue{1,N}(Tuple(c)))
end

# Trimmed :AFW are trimmed :BMM times |α|/α!
oracle_eval_afw(f::Vector{Int}, e::Tuple{Int,Int}, α::Vector{Int}, x) =
    multinomial(α...) * oracle_eval(f, e, α, x)

# ── Affine vertex-permutation map of the reference simplex ────────────────────

# Vertex V_i = the point with barycentric λⁱ = 1 (Gridap convention:
# λ = (1−Σx, x…), so V_1 = 0, V_{j+1} = e_j). Guarded by an assertion below.
ref_vertices(D) = [ [i == j+1 ? 1.0 : 0.0 for j in 1:D] for i in 1:D+1 ]

function affine_map(p::Vector{Int}, D::Int)
    V = ref_vertices(D)
    A = x -> begin
        λ = to_barycentric(x)
        y = zeros(D)
        for i in 1:D+1, j in 1:D
            y[j] += λ[i] * V[p[i]][j]
        end
        Point(y...)
    end
    x0 = A(Point(zeros(D)...))
    J  = zeros(D, D)
    for j in 1:D
        e = zeros(D); e[j] = 1.0
        J[:, j] = collect(A(Point(e...)).data) .- collect(x0.data)
    end
    A, J
end

# max |lhs − rhs| over all basis functions and points, for map permutation a:
#   lhs = (A_a^* w_μ)(x) = J_a' · w_μ(A_a x),   rhs = Σ_ν C[μ,ν] w_ν(x)
function pullback_error(b, π, a, pts, D, ev=oracle_eval)
    entries = bubble_entries(b)
    C = rotation_change_of_basis(b, π)
    A, J = affine_map(a, D)
    err = 0.0
    for (μ, t) in enumerate(entries)
        for x in pts
            lhs = transpose(J) * collect(ev(t..., A(x)).data)
            rhs = zeros(D)
            for (ν, s) in enumerate(entries)
                c = C[μ, ν]
                c == 0.0 && continue
                rhs .+= c .* collect(ev(s..., x).data)
            end
            err = max(err, maximum(abs.(lhs .- rhs)))
        end
    end
    err
end

test_points(D) = D == 2 ?
    [Point(0.13,0.27), Point(0.4,0.1), Point(0.05,0.6), Point(0.3,0.3), Point(0.45,0.45)] :
    [Point(0.1,0.2,0.3), Point(0.25,0.25,0.25), Point(0.05,0.5,0.1), Point(0.3,0.1,0.15)]

@testset "vertex/barycentric convention guard" begin
    for D in (2, 3), (i, V) in enumerate(ref_vertices(D))
        @test to_barycentric(Point(V...))[i] ≈ 1.0 atol=1e-14
    end
end

# The untrimmed rotation API is exercised through BarycentricPΛBasis 1-forms in
# both flavors; each flavor needs the oracle for its own direction form.
const ALL_BASES = (
  ((V,T,r) -> BarycentricPΛBasis(V, T,r,1; flavor=:AFW), "barycentric :AFW", oracle_eval_afw),
  ((V,T,r) -> BarycentricPΛBasis(V, T,r,1; flavor=:BMM), "barycentric :BMM", oracle_eval),
  ((V,T,r) -> BarycentricPmΛBasis(V,T,r,1; flavor=:AFW), "trimmed :AFW",     oracle_eval_afw),
  ((V,T,r) -> BarycentricPmΛBasis(V,T,r,1; flavor=:BMM), "trimmed :BMM",     oracle_eval))

@testset "direction: A_{π⁻¹} matches, A_π fails (3-cycle, D=2, r=2)" begin
    π = [2, 3, 1]
    for (make, name, ev) in ALL_BASES
        b = make(Val(2), Float64, 2)
        pts = test_points(2)
        @test pullback_error(b, π, invperm(π), pts, 2, ev) < 1e-10
        @test pullback_error(b, π, π, pts, 2, ev) > 1e-3
    end
end

@testset "C(π) == numeric pullback along A_{π⁻¹}, all π" begin
    for (make, name, ev) in ALL_BASES
        @testset "$name D=$D r=$r" for (D, rs) in ((2, (1, 2, 3)), (3, (1, 2))), r in rs
            b = make(Val(D), Float64, r)
            pts = test_points(D)
            for π in permutations(1:D+1)
                @test pullback_error(b, collect(π), invperm(collect(π)), pts, D, ev) < 1e-10
            end
        end
    end
end

@testset "RotationCache memoises and matches rotate_basis_function" begin
    for (make, name, ev) in ALL_BASES
        b  = make(Val(2), Float64, 2)
        rc = RotationCache(b)
        for π in permutations(1:3)
            m = rotation_map(rc, collect(π))
            @test m == [rotate_basis_function(b, w, collect(π)) for w in 1:length(b)]
            @test rotation_map(rc, collect(π)) === m   # cache hit
        end
        @test length(rc.maps) == 6
    end
end

end # module
