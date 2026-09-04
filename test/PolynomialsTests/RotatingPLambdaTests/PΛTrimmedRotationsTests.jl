module PΛTrimmedRotationsTests
# Change of basis on the trimmed P_r⁻Λ¹ basis induced by a vertex relabeling
# (rotation) π : ξ → λ, λ = π(ξ)
# (src/Polynomials/RotatingPLambda/PΛTrimmedRotations.jl).

using Gridap.Polynomials
using Gridap.Polynomials: rotate_multiindex, bubble_entries, rotate_basis_function
using Gridap.Polynomials: rotation_change_of_basis, trimmed_pair_sign, trimmed_pair_sort
using Combinatorics: multinomial
using LinearAlgebra
using Random
using Test

Random.seed!(1)

# ── Building-block identities ─────────────────────────────────────────────────
# Implemented directly against the formulas (ξ^β, ϕ), not against the rotation
# code itself.

# ϕ(λ;e1,e2) = λ^{e1} dλ^{e2} − λ^{e2} dλ^{e1}, as a coefficient vector in the
# ambient frame dλ¹,…,dλᴺ.
function ϕ(e1::Int, e2::Int, N::Int, λ)
  c = zeros(eltype(λ), N)
  c[e2] += λ[e1]
  c[e1] -= λ[e2]
  c
end

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

@testset "ϕ pullback: π⁻¹*(ϕ(ξ;e1,e2)) = ϕ(λ;π(e1),π(e2))" begin
    N = 3
    pairs = [(1,2),(2,1),(1,3),(3,1),(2,3),(3,2)]
    for trial in 1:10, (e1,e2) in pairs
        π   = randperm(N)
        ξ   = rand(N)
        λ   = ξ[invperm(π)]
        c_pb   = ϕ(e1, e2, N, ξ)[invperm(π)]
        c_dir  = ϕ(π[e1], π[e2], N, λ)
        @test c_pb ≈ c_dir
    end
end

@testset "ε/↑ split: π⁻¹*(ϕ(ξ;e)) = ε(π(e))*ϕ(λ;π(e)↑)" begin
    N = 3
    pairs = [(1,2),(2,1),(1,3),(3,1),(2,3),(3,2)]
    for trial in 1:10, (e1,e2) in pairs
        π   = randperm(N)
        ξ   = rand(N)
        λ   = ξ[invperm(π)]
        c_pb = ϕ(e1, e2, N, ξ)[invperm(π)]
        πe1,πe2 = π[e1], π[e2]
        ε = trimmed_pair_sign(πe1,πe2)
        eπ1,eπ2 = trimmed_pair_sort(πe1,πe2)
        c_split = ε .* ϕ(eπ1, eπ2, N, λ)
        @test c_pb ≈ c_split
    end
end

@testset "Hit identity: ξ_v0·ϕ(v1,v2) = ξ_v1·ϕ(v0,v2) - ξ_v2·ϕ(v0,v1)" begin
    N = 4
    triples = [(1,2,3),(2,1,3),(1,3,4),(4,2,1),(2,3,4)]
    for trial in 1:10, (v0,v1,v2) in triples
        ξ = rand(N)
        lhs = ξ[v0] .* ϕ(v1,v2,N,ξ)
        rhs = ξ[v1] .* ϕ(v0,v2,N,ξ) .-
              ξ[v2] .* ϕ(v0,v1,N,ξ)
        @test lhs ≈ rhs
    end
end

# ── Basis-level ───────────────────────────────────────────────────────────────

trimmed_basis(D, r; flavor=:BMM) = BarycentricPmΛBasis(Val(D), Float64, r, 1; flavor)

@testset "trimmed basis dimension sanity (D=2)" begin
    @test length(trimmed_basis(2, 1)) == 3  # lowest-order Whitney edges
    @test length(trimmed_basis(2, 2)) == 8
end

D, r = 2, 2
n = length(trimmed_basis(D, r))

@testset "flavor=$flavor" for flavor in (:BMM, :AFW)
    b = trimmed_basis(D, r; flavor)

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

# The two flavors scale basis function w by multinomial(α_w), so their change-of-
# basis matrices are conjugate by that diagonal. This pins the hit weights: a
# round trip alone cannot, since S·C·S⁻¹ round trips for ANY diagonal S.
@testset "flavors are diagonally conjugate: C_AFW = S C_BMM S⁻¹" begin
    for (D, r) in ((2,2), (2,3), (3,2), (3,3))
        bmm = trimmed_basis(D, r; flavor=:BMM)
        afw = trimmed_basis(D, r; flavor=:AFW)
        @test bubble_entries(afw) == bubble_entries(bmm)

        S = Diagonal([float(multinomial(α...)) for (_, _, α) in bubble_entries(bmm)])
        swap, rev, cycle = [2, 1, 3:D+1...], collect(D+1:-1:1), [2:D+1..., 1]
        for π in (swap, rev, cycle)
            @test rotation_change_of_basis(afw, π) ≈
                  S * rotation_change_of_basis(bmm, π) * inv(S)
        end
    end
end

end # module
