module BarycentricPΛBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Helpers
using Combinatorics: combinations, multiexponents, multinomial, permutations
using ForwardDiff
using LinearAlgebra
using Random
using StaticArrays

using Gridap.Fields: Point, return_cache, evaluate!
using Gridap.Polynomials: bernstein_term_id
using Gridap.Polynomials: bubble_entries, rotate_basis_function, rotate_face_set
using Gridap.Polynomials: rotate_multiindex, rotation_change_of_basis
using Gridap.Polynomials: trimmed_pair_sign, trimmed_pair_sort
using Gridap.Polynomials: RotationCache, rotation_map

r = 3 # all possible bubble spaces are non empty

# Combination ordering validation
for D in 1:6
  for k in 0:D
    # standard lexicographic order (left-to-right)
    acc = true
    for (I_ind, I) in enumerate(sorted_combinations(D, k))
      acc &= combination_index(I, D) == I_ind
    end
    @test acc || "A combination index was wrong for $D, $k"
    # right-to-left lexicographic order
    right_to_left = true
    acc = true
    for (I_ind, I) in enumerate(sorted_combinations(D, k; right_to_left))
      acc &= combination_index(I, D; right_to_left) == I_ind
    end
    @test acc || "A combination index was wrong for $D, $k, right-to-left order"
  end
end


# Bubble indices validation

D = 3
N = D+1
for k in 0:D
  w_prev = 0
  for (F, bubble_functions) in PmΛ_bubbles(r,k,D)
    @test issorted(F)  &&  (k ≤ length(F) ≤ N)  &&  (F ⊆ 1:N) || (k,D,F)

    passed = true
    for (w, α, α_id, J, sub_J_ids, sup_α_ids) in bubble_functions
      passed == passed && w == w_prev+1
      w_prev = w

      passed = passed && all(α .≥ 0) && sum(α)==r-1 && length(α)==N &&
            α_id == bernstein_term_id(α) &&
            all( bernstein_term_id( [α[j]+Int(i==j) for j in eachindex(α)] ) == αpi_id
                for (i,αpi_id) in enumerate(sup_α_ids) )

      passed = passed && issorted(J) && length(J)==k+1 && (J ⊆ 1:N) &&
            all( combination_index(J[J .≠ J[i]], N) == Jsi_id for (i,Jsi_id) in enumerate(sub_J_ids) )
    end
    @test passed || (r, k, D, F, bubble_functions)
  end
  @test w_prev == binomial(r+k-1,k)*binomial(D+r,D-k)
end

for k in 0:D
  w_prev = 0
  for (F, bubble_functions) in PΛ_bubbles(r,k,D)
    @test issorted(F)  &&  (k ≤ length(F) ≤ N)  &&  (F ⊆ 1:N) || (k,D,F)

    passed = true
    for (w, α, α_id, J) in bubble_functions
      passed == passed && w == w_prev+1
      w_prev = w

      passed = passed && all(α .≥ 0) && sum(α)==r && length(α)==N &&
            α_id == bernstein_term_id(α)

      passed = passed && issorted(J) && length(J)==k && (J ⊆ 1:N)
    end
    @test passed || (r, k, D, F, bubble_functions)
  end
  @test w_prev == binomial(r+k,k)*binomial(D+r,D-k)
end


# Bases tests

using Gridap.Polynomials: _minusone_if_even_else_one, _findfirst_val_or_zero

function _test_testvalue(b, Bx, Gx, Hx)
  b0 = testvalue(b)
  @test b0 isa typeof(b)

  @test evaluate(b0,x)                   isa typeof(Bx)
  @test evaluate(Broadcasting(∇)(b0),x)  isa typeof(Gx)
  @test evaluate(Broadcasting(∇∇)(b0),x) isa typeof(Hx)
end

#############################################################################
# Former logic with analytical computation of the basis in the Reference FE #
#############################################################################

function _test_reference_basis(b,D,r,k,flavor)
  flavor != :AFW && return
  V = value_type(b)
  V <: Real && (V = VectorValue{1,V})

  if b isa BarycentricPmΛBasis
    LN = binomial(D+1,k)
    m = zeros(V, LN)
    _compute_PmΛ_basis_reference_coefficients!(m,k,D,b._indices)
    @test all(@. norm(b.m - m) < 1.e-15)

  else       #BarycentricPΛBasis
    Ψ = zeros(V, length(b))
    _compute_PΛ_basis_reference_form_coefficient!(Ψ,r,k,b._indices)
    @test all(@. norm(b.Ψ - Ψ) < 1.e-15)
  end
end

function _compute_PmΛ_basis_reference_coefficients!(m,k,D,indices)

  if iszero(k) # so V is scalar, no change of basis
    m .= 1
    return nothing
  end

  V = eltype(m)
  m_J = Mutable(V)(undef)
  @inbounds for (J_id, J) in enumerate(sorted_combinations(D+1,k))
    s = Int(isone(J[1]))
    for (I_id, I, I_sgn) in indices.components
      n = count(i-> (J[i]-1)∉I, (1+s):k)
      if iszero(n)
        p = _findfirst_val_or_zero(j-> (I[j]+1)∉J, 1, k)
        m_J[I_id] = I_sgn*_minusone_if_even_else_one(p+1)
      else
        m_J[I_id] = 0
      end
    end
    m[J_id] = m_J
  end
  nothing
end

function _compute_PΛ_basis_reference_form_coefficient!(Ψ,r,k,indices)

  """
      _hat_Ψ(r,::Val{k},α,F,I,J,T)::T

  BarycentricPΛBasis.Ψ matrix elements in the reference simplex, T is the scalar return type

  This is actually not faster than computing the matrices and the minors
  explicitely like when vertices are given, but might be usefull in case we want
  to compute these at compile time one day.
  """
  function _hat_Ψ(r,Vk::Val{k},α,F,I,J,::Type{T})::T where {T,k}
    @check sum(α) == r
    @check length(I) == length(J) == k

    iszero(k) && return one(T) # 0 forms

    _u(i::Int,F,I)   = Int(isone(F[1])) - Int(I[i]+1 in F)
    _u(F::Vector{Int},I,Vk) = ntuple(i->_u(i,F,I), Vk)
    _v(j::Int,α,J,r) = α[J[j]]/r
    _v(α::Vector{Int},J,r,Vk) = ntuple(j->_v(j,α,J,r), Vk)

    @inbounds begin

      s = Int(isone(J[1]))
      n = count(i-> (J[i]-1)∉I, (1+s):k)

      n > 1 && return 0. # rank M_IJ inferior to 2

      p = _findfirst_val_or_zero(j-> (I[j]+1)∉J, 1, k)

      if isone(n)        # rank M_IJ is 1
        m = _findfirst_val_or_zero(i-> (J[i]-1)∉I, (s+1), k)
        u_p, v_m = _u(p,F,I), _v(m,α,J,r)
        sgn = _minusone_if_even_else_one(m+p+1)
        iszero(s) && return sgn*u_p*v_m

        q = _findfirst_val_or_zero(j-> (I[j]+s)∉J, (p+1), k)
        u_q = _u(q,F,I)
        sgn *= _minusone_if_even_else_one(q+1)
        return sgn * v_m * (u_q - u_p)
      end

      u, v = _u(F,I,Vk), _v(α,J,r,Vk)
      if iszero(s)
        return 1 + sum( u .* v )
      else
        Ψ_IJ = one(T)
        sum_v = v[1]
        for l in 1:p-1
          vlp = v[l+1]
          sum_v += vlp
          Ψ_IJ += vlp*u[l]
        end
        for l in (p+1):k
          vl = v[l]
          sum_v += vl
          Ψ_IJ += vl*u[l]
        end
        sgn = _minusone_if_even_else_one(p+1)
        return sgn * (Ψ_IJ - u[p]*sum_v)
      end

    end
    @unreachable
  end

  Vk = Val(k)
  V = eltype(Ψ)
  T = eltype(V)
  Ψw = Mutable(V)(undef)

  iszero(r) && return _order_0_Ψ!(Ψ)

  @inbounds for (F, bubble_functions) in indices.bubbles
    for (w, α, _, J) in bubble_functions
      for (I_id, I, I_sgn) in indices.components
        Ψw[I_id] = I_sgn * _hat_Ψ(r,Vk,α,F,I,J,T)
      end
      Ψ[w] = Ψw
    end
  end
  nothing
end

##############################
# testing function and tests #
##############################

function _test_basis(VD::Val{D}, T, r, k, vertices) where D
  for PΛB in (BarycentricPmΛBasis, BarycentricPΛBasis)
    for flavor in (:AFW, :BMM)
      b   = PΛB(VD,T,r,k;flavor)
      @test contains(sprint(show, MIME"text/plain"(), b._indices), "Λᵏ(△ᴰ) basis indices, r=$r k=$k D=$D")
      @test_nowarn print_indices(b,IOBuffer())
      @test get_orders(b) == tfill(r,Val(D))

      b2  = PΛB(VD,T,r,k; indices=b._indices, flavor) # indices recycling
      @test b == b2
      @test b2._indices == b._indices

      faces = [bubble[1] for bubble in get_bubbles(b)] # bubble space selection
      b2  = PΛB(b, faces...)
      @test b == b2

      _test_reference_basis(b,D,r,k,flavor)

      Bx = evaluate(b,x)
      Gx = evaluate(Broadcasting(∇)(b),x)
      Hx = evaluate(Broadcasting(∇∇)(b),x)
      _test_testvalue(b, Bx, Gx, Hx)

      bv  = PΛB(VD,T,r,k,vertices; flavor)
      Bx = evaluate(bv,x)
      Gx = evaluate(Broadcasting(∇)(bv),x)
      Hx = evaluate(Broadcasting(∇∇)(bv),x)
      _test_testvalue(bv, Bx, Gx, Hx)
    end
  end
end

T = Float64

# 0D                                           0D #
D = 0
vertices = [Point{D,T}()]
x = [vertices[1]]
k = 0

_test_basis(Val(D), T, r, k, vertices)

# 1D                                           1D #
D = 1
Pt = Point{D,T}
vertices = [Pt(0.),Pt(1.)]
x = [xi for xi in vertices]

for k in 0:D
  _test_basis(Val(D), T, r, k, vertices)
end


# 2D                                           2D #
D = 2
Pt = Point{D,T}
vertices = [Pt(0., 0),Pt(1.,0),Pt(0.,1.)]
x = [xi for xi in vertices]

for k in 0:D
  _test_basis(Val(D), T, r, k, vertices)
end

# 3D                                           3D #
D = 3
Pt = Point{D,T}
vertices = [Pt(0.,0,0),Pt(1.,0,0),Pt(0,1.,0),Pt(0,0,1.)]
x = [xi for xi in vertices]

for k in 0:D
  _test_basis(Val(D), T, r, k, vertices)
end

# 4D                                           4D #
D = 4
Pt = Point{D,T}
vertices = [Pt(0.,0,0,0),Pt(1.,0,0,0),Pt(0,1.,0,0),Pt(0,0,1.,0),Pt(0,0,0,1.)]
x = [xi for xi in vertices]

for k in 0:D
  k == 2 && continue # no vector proxy
  _test_basis(Val(D), T, r, k, vertices)
end


@testset "PΛRotations" begin
  # Change of basis on BarycentricPΛBasis 1-forms induced by a vertex relabeling
  # (rotation) π : ξ → λ, λ = π(ξ)  (src/Polynomials/RotatingPLambda/PΛRotations.jl).

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

  FLAVORS = (:AFW, :BMM)

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
end

@testset "PΛTrimmedRotations" begin
  # Change of basis on the trimmed P_r⁻Λ¹ basis induced by a vertex relabeling
  # (rotation) π : ξ → λ, λ = π(ξ)
  # (src/Polynomials/RotatingPLambda/PΛTrimmedRotations.jl).

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
end

@testset "PΛRotationPullback" begin
  # End-to-end validation of the rotation change of basis against a NUMERIC
  # pullback, for the untrimmed BarycentricPΛBasis 1-forms and the trimmed
  # BarycentricPmΛBasis ones, in both flavors each.
  #
  # The self-consistency tests in the PΛRotations and PΛTrimmedRotations testsets
  # check the closed-form index calculus against itself; this one closes the
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

  # ── Inlined helpers (from the ExteriorGridap.jl SimplexFEECBases oracle) ──────

  # Gridap convention: λ = (1−Σx, x…)
  to_barycentric(x::Point{D,T}) where {D,T} = (one(T) - sum(x.data), x.data...)

  # Project an ambient barycentric 1-form (dλ¹,…,dλᴺ) to the physical D = N−1
  # frame by imposing dλ¹ = −dλ² − … − dλᴺ: phys[j] = c[j+1] − c[1].
  reduce_ambient(ω::ExteriorFormValue{1,N}) where N =
      ExteriorFormValue{1,N-1}(ntuple(j -> ω.data[j+1] - ω.data[1], N-1))

  # ── Oracles (identical formulas to the basis testsets) ───────────────────────

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
  ALL_BASES = (
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
end

@testset "PΛBasisGradients" begin
  # Gradient path of the BarycentricPΛBasis and BarycentricPmΛBasis 1-forms (∇
  # via Gridap's FieldGradientArray machinery) validated against central finite
  # differences of the value path.
  # Convention: ∇u[a,j] = ∂u_j/∂x_a (TensorValue{D,D}).

  test_points(D) = D == 2 ?
      [Point(0.13,0.27), Point(0.3,0.3), Point(0.45,0.45), Point(0.1,0.1)] :
      [Point(0.1,0.2,0.3), Point(0.25,0.25,0.25), Point(0.3,0.1,0.15)]

  perturb(x, a, h, D) = Point(ntuple(j -> x[j] + (j == a ? h : 0.0), D)...)

  @testset "∇(basis) == finite differences" begin
      h = 1e-5
      for (make, name) in (((V,T,r) -> BarycentricPΛBasis(V,T,r,1;  flavor=:AFW), "untrimmed :AFW"),
                           ((V,T,r) -> BarycentricPΛBasis(V,T,r,1;  flavor=:BMM), "untrimmed :BMM"),
                           ((V,T,r) -> BarycentricPmΛBasis(V,T,r,1; flavor=:AFW), "trimmed :AFW"),
                           ((V,T,r) -> BarycentricPmΛBasis(V,T,r,1; flavor=:BMM), "trimmed :BMM"))
          @testset "$name D=$D r=$r" for (D, rs) in ((2, (1, 2, 3)), (3, (1, 2))), r in rs
              b   = make(Val(D), Float64, r)
              pts = test_points(D)
              vals  = evaluate(b, pts)
              grads = evaluate(Broadcasting(∇)(b), pts)
              for (i, x) in enumerate(pts), w in 1:length(b)
                  G = grads[i, w]
                  for a in 1:D
                      vp = evaluate(b, [perturb(x, a,  h, D)])[1, w]
                      vm = evaluate(b, [perturb(x, a, -h, D)])[1, w]
                      for j in 1:D
                          fd = (vp.data[j] - vm.data[j]) / (2h)
                          @test isapprox(G[a, j], fd; atol=1e-6)
                      end
                  end
              end
          end
      end
  end
end

@testset "RotatingPΛBases" begin
  # The rotating P_rΛ¹ basis, i.e. BarycentricPΛBasis 1-forms with flavor=:BMM,
  # built from the directional 1-form
  #
  #   ψ(f,k,s(α)) = dλᵏ − (𝟙[k∈supp(α)]/|supp(α)|) Σ_{i∈f} dλⁱ
  #
  # A self-contained, hand-rolled oracle (spanning set + filter + numeric
  # evaluation, independent of the package internals) is implemented below and
  # cross-checked against BarycentricPΛBasis built with flavor=:BMM.

  # ── Inlined oracle helpers ────────────────────────────────────────────────────

  # Gridap convention: λ = (1−Σx, x…)
  to_barycentric(x::Point{D,T}) where {D,T} = (one(T) - sum(x.data), x.data...)

  # Project an ambient barycentric 1-form (dλ¹,…,dλᴺ) to the physical D = N−1
  # frame by imposing dλ¹ = −dλ² − … − dλᴺ: phys[j] = c[j+1] − c[1].
  reduce_ambient(ω::ExteriorFormValue{1,N}) where N =
      ExteriorFormValue{1,N-1}(ntuple(j -> ω.data[j+1] - ω.data[1], N-1))

  # ψ(f,k,α) = dλᵏ − (𝟙[k∈supp(α)] / |supp(α)|) · Σ_{i∈f} dλⁱ
  # returned as a coefficient vector in the AMBIENT (D+1)-dim barycentric frame
  # dλ¹,…,dλ^{D+1}  (reduce_ambient later projects this to the D physical dλ's).
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
  #   degree-≤r polynomials via Σλⁱ=1),  f a face of Δ_D with dim(f) ≥ K,  k ∈ f,
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
      ω   = ExteriorFormValue{1,N}(Tuple(rotating_psi(f, k, α)))
      val * reduce_ambient(ω)
  end

  # ── Driver: D=2, K=1, r=3 ──────────────────────────────────────────────────────
  # r must be ≥ 3 for ψ to differ from the AFW direction form: below that, α_j/|α|
  # and s(α)_j/|s(α)| agree on every basis function.
  D, K, r = 2, 1, 3
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
  # BarycentricPΛBasis{D} with flavor=:BMM generates the same (f,k,α) triples (the
  # vertex k of the oracle is the single index J[1] of a bubble function) and the
  # same ψ formula, but contracted against the simplex's actual ∂λ/∂x Jacobian
  # instead of going through reduce_ambient. On the default (unit) simplex the two
  # contractions coincide, so the coefficients must agree exactly.
  #
  # The basis is vector proxied: for k=1 the proxy is the identity, so its Ψ holds
  # the same dx components as the oracle's reduced 1-form.
  bmm_basis(r) = BarycentricPΛBasis(Val(D), Float64, r, K; flavor=:BMM)

  # (f,k,α) → w, the oracle's key onto the package's basis function index
  function bmm_index(b)
      index = Dict{Tuple{Vector{Int},Int,Vector{Int}},Int}()
      for (F, bfs) in get_bubbles(b), (w, α, _, J) in bfs
          index[(F, J[1], α)] = w
      end
      index
  end

  @testset "BarycentricPΛBasis(:BMM) == hand-rolled oracle coefficients" begin
      pkg_basis = bmm_basis(r)
      @test length(pkg_basis) == length(basis)

      pkg_index = bmm_index(pkg_basis)

      @testset "(f,k,α) triples present and ψ values match" for (f, k, α) in basis
          key = (f, k, collect(α))
          @test haskey(pkg_index, key)
          pkg_ψ    = collect(pkg_basis.Ψ[pkg_index[key]].data)
          oracle_ψ = collect(reduce_ambient(ExteriorFormValue{1,D+1}(Tuple(rotating_psi(f, k, α)))).data)
          @test pkg_ψ ≈ oracle_ψ atol=1e-12
      end
  end

  # ── Evaluation path: return_cache / evaluate! of the package basis ───────────
  #
  # The coefficient cross-check above never exercises the Gridap evaluation
  # machinery (_return_cache/_setsize!/_evaluate_nd!). Evaluate the package basis
  # pointwise and compare against the oracle.
  @testset "BarycentricPΛBasis(:BMM) evaluate! == oracle, pointwise" begin
      pkg_basis = bmm_basis(r)
      pkg_index = bmm_index(pkg_basis)
      cache    = return_cache(pkg_basis, pts)
      pkg_vals = evaluate!(cache, pkg_basis, pts)   # (np, ndof) matrix
      for (f, k, α) in basis, (i, p) in enumerate(pts)
          w = pkg_index[(f, k, collect(α))]
          @test collect(pkg_vals[i, w].data) ≈ collect(rotating_eval(f, k, α, p).data) atol=1e-12
      end
  end
end

@testset "TrimmedPΛBases" begin
  # The trimmed P_r⁻Λ¹ basis, built from the Whitney (directional) 1-form
  #
  #   ϕ(ξ;e1,e2) = ξ^{e1} dξ^{e2} − ξ^{e2} dξ^{e1}
  #
  # A self-contained, hand-rolled oracle (spanning set + filter + numeric
  # evaluation, independent of the package internals) is implemented below and
  # cross-checked against BarycentricPmΛBasis with flavor=:BMM. Unlike ψ in the
  # untrimmed case (the RotatingPΛBases testset), ϕ is NOT constant-coefficient —
  # it is linear in λ — so the cross-check compares pointwise numeric evaluations
  # rather than constant coefficient vectors.

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
end

end # module
