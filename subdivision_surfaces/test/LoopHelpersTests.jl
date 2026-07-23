module LoopHelpersTests

using Test
using LinearAlgebra
using GridapSubdivisionSurfaces.ReferenceFEs.LoopHelpers

@testset "_jordan_fix_N3 hardcoded values match the general derivation" begin
  # independent re-derivation of the N=3 Jordan-chain fix (see
  # `_jordan_fix_N3` in LoopHelpers.jl), used only to validate the hardcoded
  # values there haven't drifted from the general algorithm they came from
  u0 = U0(3)[:,2]
  λ = Δ()[4,4]
  s12 = S12()
  rhs = S11(3)*u0
  Mj = λ*I - s12
  leftnull = nullspace(Matrix(Mj'))
  W16 = W1()[:,4:5]
  ab = (leftnull'*W16) \ (leftnull'*rhs)
  wa = W16*ab
  perp = ab/norm(ab)
  wb = W16*[-perp[2], perp[1]]
  u1 = pinv(Mj)*(rhs - wa)

  u1_hc, wa_hc, wb_hc = LoopHelpers._jordan_fix_N3()
  @test u1 ≈ u1_hc atol=1e-14
  @test wa ≈ wa_hc atol=1e-14
  # wb is only defined up to scale (any nonzero vector spanning the
  # complementary eigenspace direction is a valid eigenvector), so check
  # proportionality rather than equality
  @test rank(hcat(wb, wb_hc)) == 1
end

@testset "A = V*Λ*V_inv, N=3:20" begin
  for N in 3:20
    A_N, V_N, V_inv_N, Λ_N = A(N), V(N), V_inv(N), Λ(N)
    @test norm(A_N*V_N - V_N*Λ_N) < 1e-14
    @test norm(V_N*V_inv_N - I) < 1e-14
  end
end

@testset "P(N,k) and A_bar(N) consistency" begin
  # A_bar's first K rows are exactly A, by construction
  for N in 3:10
    K = N+6
    @test A_bar(N)[1:K, :] == A(N)
  end

  # for N >= 4, each P(N,k) picks 12 distinct, in-range columns, each with
  # weight 1, and P(N,k)*A_bar(N) is a partition of unity (rows sum to 1)
  for N in 4:10, k in 1:3
    Pk = P(N, k)
    @test size(Pk) == (12, N+12)
    cols = [findfirst(!=(0.0), Pk[r,:]) for r in 1:12]
    @test allunique(cols)
    @test all(c -> 1 <= c <= N+12, cols)
    @test all(r -> Pk[r, cols[r]] == 1.0, 1:12)

    B = Pk * A_bar(N)
    @test all(isapprox.(sum(B, dims=2), 1.0))
  end

  # N=3 is a documented exception: P(3,2) picks the same physical vertex at
  # two different canonical positions (see P's docstring); P(3,1)/P(3,3) are
  # unaffected. The partition-of-unity property still holds regardless.
  for k in (1, 3)
    cols = [findfirst(!=(0.0), P(3,k)[r,:]) for r in 1:12]
    @test allunique(cols)
  end
  cols2 = [findfirst(!=(0.0), P(3,2)[r,:]) for r in 1:12]
  @test !allunique(cols2)
  @test length(unique(cols2)) == 11

  for k in 1:3
    B3 = P(3,k) * A_bar(3)
    @test all(isapprox.(sum(B3, dims=2), 1.0))
  end
end

end # module
