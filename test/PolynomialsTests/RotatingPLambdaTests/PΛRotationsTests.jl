module PΛRotationsTests
# Change of basis on BarycentricPΛBasis 1-forms induced by a vertex relabeling
# (rotation) π : ξ → λ, λ = π(ξ)  (src/Polynomials/RotatingPLambda/PΛRotations.jl).

using Gridap.Polynomials
using Gridap.Polynomials: rotate_multiindex, rotate_face_set, bubble_entries
using Gridap.Polynomials: rotate_basis_function, rotation_change_of_basis
using LinearAlgebra
using Random
using Test

Random.seed!(1)

# The direction 1-form φ(F,k,α) of the basis, on the ambient dλ¹,…,dλᴺ coframe:
#   dλᵏ - (αₖ/|α|) Σ_{i∈F} dλⁱ        for :AFW
#   dλᵏ - (s(α)ₖ/|s(α)|) Σ_{i∈F} dλⁱ  for :BMM, s(α) the support indicator
# Written out here rather than taken from the package, so the identities below
# check the rotation code against the formulas instead of against themselves.
function ambient_phi(F::Vector{Int}, k::Int, α::Vector{Int}, N::Int, flavor::Symbol)
    s = flavor === :BMM ? [Int(a > 0) for a in α] : α
    ρ = sum(s)
    c = zeros(N)
    c[k] += 1.0
    iszero(ρ) || for i in F; c[i] -= s[k]/ρ; end
    c
end

const FLAVORS = (:AFW, :BMM)

# ── Building-block identities ─────────────────────────────────────────────────
# Implemented directly against the formulas (ξ^β, ambient_phi), not against the
# rotation code itself, so these are a genuine check of that file rather than a
# tautology.

@testset "BB pullback: π⁻¹*(ξ^β) = λ^π(β)" begin
    N = 3
    for trial in 1:10
        π  = randperm(N)
        β  = rand(0:3, N)
        ξ  = rand(N)
        λ  = ξ[invperm(π)]   # λⁱ = ξ^{π⁻¹(i)}
        ξβ  = prod(ξ[i]^β[i] for i in 1:N)
        πβ  = rotate_multiindex(β, π)
        λπβ = prod(λ[i]^πβ[i] for i in 1:N)
        @test ξβ ≈ λπβ
    end
end

@testset "φ pullback: π⁻¹*(φ(ξ;F,k,α)) = φ(λ;π(F),π(k),π(α)) [$flavor]" for flavor in FLAVORS
    N = 3
    triples = [([1,2],1,[0,2,0]), ([1,2],2,[1,1,0]), ([1,3],3,[2,0,0]),
               ([1,2,3],1,[0,1,1]), ([1,2,3],3,[1,1,0])]
    for trial in 1:10, (F,k,α) in triples
        π    = randperm(N)
        c    = ambient_phi(F, k, α, N, flavor)
        cpb  = c[invperm(π)]   # pullback of dξ-coefficients to dλ
        πF, πk, πα = rotate_face_set(F, π), π[k], rotate_multiindex(α, π)
        cdir = ambient_phi(πF, πk, πα, N, flavor)
        @test cpb ≈ cdir
    end
end

@testset "Filter-hit re-expansion: φ(F,min(F),α) = -∑_{vi∈F∖min(F)} φ(F,vi,α) [$flavor]" for flavor in FLAVORS
    N = 3
    full_support_triples = [([1,2],[1,1,0]), ([1,3],[1,0,1]), ([2,3],[0,1,1]),
                             ([1,2,3],[1,1,1])]
    for (F,α) in full_support_triples
        k0  = minimum(F)
        lhs = ambient_phi(F, k0, α, N, flavor)
        rhs = -sum(ambient_phi(F, vi, α, N, flavor) for vi in F if vi != k0)
        @test lhs ≈ rhs
    end
end

# ── Basis-level change of basis ───────────────────────────────────────────────

D, r = 2, 2

@testset "basis-level change of basis [$flavor]" for flavor in FLAVORS
    b = BarycentricPΛBasis(Val(D), Float64, r, 1; flavor=flavor)
    n = length(b)

    @testset "rotation_change_of_basis: identity permutation gives I" begin
        C = rotation_change_of_basis(b, collect(1:D+1))
        @test C ≈ Matrix(1.0I, n, n)
    end

    @testset "rotation_change_of_basis: round trip via invperm(π), no matrix inverse" begin
        perms = [[2,1,3], [1,3,2], [3,2,1], [2,3,1], [3,1,2]]
        for π in perms
            C    = rotation_change_of_basis(b, π)
            Cinv = rotation_change_of_basis(b, invperm(π))
            @test C * Cinv ≈ Matrix(1.0I, n, n)
            @test Cinv * C ≈ Matrix(1.0I, n, n)
        end
    end

    @testset "rotation_change_of_basis rows agree with rotate_basis_function" begin
        π = [2,3,1]
        C = rotation_change_of_basis(b, π)
        for w in 1:n
            row = zeros(n)
            for (c, w′) in rotate_basis_function(b, w, π)
                row[w′] += c
            end
            @test row ≈ C[w, :]
        end
    end
end

@testset "the rotation API rejects k ≠ 1" begin
    @test_throws Exception bubble_entries(BarycentricPΛBasis(Val(2), Float64, 2, 2))
end

end # module
