module PΛRotationsTests
# Change of basis on RotatingPΛBasis induced by a vertex relabeling (rotation)
# π : ξ → λ, λ = π(ξ)  (src/Polynomials/RotatingPLambda/PΛRotations.jl).

using Gridap.Polynomials
using Gridap.Polynomials: _rotating_ambient_psi
using LinearAlgebra
using Random
using Test

Random.seed!(1)

# ── Building-block identities ─────────────────────────────────────────────────
# Implemented directly against the formulas (ξ^β, _rotating_ambient_psi), not
# against the rotation code itself, so these are a genuine check of that file
# rather than a tautology.

@testset "BB pullback: π⁻¹*(ξ^β) = λ^π(β)" begin
    N = 3
    for trial in 1:10
        π  = randperm(N)
        β  = rand(0:3, N)
        ξ  = rand(N)
        λ  = ξ[invperm(π)]   # λ_i = ξ_{π⁻¹(i)}
        ξβ  = prod(ξ[i]^β[i] for i in 1:N)
        πβ  = rotate_multiindex(β, π)
        λπβ = prod(λ[i]^πβ[i] for i in 1:N)
        @test ξβ ≈ λπβ
    end
end

@testset "ψ pullback: π⁻¹*(ψ(ξ;F,k,α)) = ψ(λ;π(F),π(k),π(α))" begin
    N = 3
    triples = [([1,2],1,[0,2,0]), ([1,2],2,[1,1,0]), ([1,3],3,[2,0,0]),
               ([1,2,3],1,[0,1,1]), ([1,2,3],3,[1,1,0])]
    for trial in 1:10, (F,k,α) in triples
        π    = randperm(N)
        c    = _rotating_ambient_psi(F, k, α, N)
        cpb  = c[invperm(π)]   # pullback of dξ-coefficients to dλ
        πF, πk, πα = rotate_face_set(F, π), π[k], rotate_multiindex(α, π)
        cdir = _rotating_ambient_psi(πF, πk, πα, N)
        @test cpb ≈ cdir
    end
end

@testset "Filter-hit re-expansion: ψ(F,min(F),α) = -∑_{vi∈F∖min(F)} ψ(F,vi,α)" begin
    N = 3
    full_support_triples = [([1,2],[1,1,0]), ([1,3],[1,0,1]), ([2,3],[0,1,1]),
                             ([1,2,3],[1,1,1])]
    for (F,α) in full_support_triples
        k0  = minimum(F)
        lhs = _rotating_ambient_psi(F, k0, α, N)
        rhs = -sum(_rotating_ambient_psi(F, vi, α, N) for vi in F if vi != k0)
        @test lhs ≈ rhs
    end
end

# ── Basis-level change of basis ───────────────────────────────────────────────

D, r = 2, 2
b = RotatingPΛBasis(Val(D), Float64, r)
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

end # module
