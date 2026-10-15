module PullbackKoszulTests
# Value-level pullback/pushforward under a linear map (Jacobian), the Koszul
# operator κ_x(ω) = ι_x(ω).

using Gridap.TensorValues
using Test
using LinearAlgebra: pinv
using Gridap.Arrays: get_array

# ── Pullback (value level) ───────────────────────────────────────────────────

# The pullback takes the transposed Jacobian Jt[j,i] = ∂φⁱ/∂ξʲ (a.k.a.
# gradient), whose chart index j is covariant (it contract with vector
# components) and ambient index i contravariant (the resulting ambiant vector
# component).

# Identity map on ℝ²: Jt = I₂, pullback should be identity
Jt_id = TensorValue{cov2,con2,Float64,4}(1.0,0.0,0.0,1.0)
ω1   = CoVectorValue{2}((3.0, 5.0))
@test pullback(ω1, Jt_id).data == (3.0, 5.0)

vol2 = ExteriorFormValue{2,2}((1.0,))
@test pullback(vol2, Jt_id).data == (1.0,)

# Scaling map ℝ² → ℝ²: φ(ξ) = (aξ¹, bξ²), Jt = diag(a,b)
a, b = 2.0, 3.0
Jt_sc = TensorValue{cov2,con2,Float64,4}(a,0.0,0.0,b)

# pullback of dx¹: (φ*dx¹)_j = Jt[j,1] → a dx¹  (only j=1 gives non-zero)
@test collect(pullback(CoVectorValue{2}((1.0,0.0)), Jt_sc).data) ≈ [a, 0.0]
@test collect(pullback(CoVectorValue{2}((0.0,1.0)), Jt_sc).data) ≈ [0.0, b]

# pullback of dx¹∧dx² under Jt_sc = diag(a,b): (φ*vol)_12 = det(Jt) = a*b
@test collect(pullback(vol2, Jt_sc).data) ≈ [a*b]   atol=1e-12

# Embedding ℝ² → ℝ³: project onto (x¹,x²,0) plane
# Jt = [1 0 0; 0 1 0], 2 chart rows and 3 ambient columns
Jt_emb = TensorValue{cov2,con3,Float64,6}(1.0,0.0,  0.0,1.0,  0.0,0.0)
ω1_3d = CoVectorValue{3}((1.0, 0.0, 0.0))  # dx¹ in ℝ³
pb    = pullback(ω1_3d, Jt_emb)
@test typeof(pb) <: CoVectorValue{cov2}
@test pb.data == (1.0, 0.0)   # dx¹ pulls back to dx¹ in ℝ²

ω1_3d_2 = CoVectorValue{3}((0.0, 1.0, 0.0))  # dx² in ℝ³
@test pullback(ω1_3d_2, Jt_emb).data == (0.0, 1.0)       # dx² in ℝ²

ω1_3d_3 = CoVectorValue{3}((0.0, 0.0, 1.0))  # dx³ in ℝ³ → 0 in ℝ²
@test collect(pullback(ω1_3d_3, Jt_emb).data) ≈ [0.0, 0.0]

# The ambient index of Jt is the contravariant one, so it is the index that
# contracts with the form. Swapping the two rejects the Jacobian rather than
# silently pulling back along the transposed map.
@test_throws ArgumentError pullback(ω1_3d, transpose(Jt_emb))
@test_throws ArgumentError pullback(ω1, TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0))
@test_throws DimensionMismatch pullback(ω1, TensorValue{cov2,con3,Float64,6}(1.,0,0,1,0,0))

# A non-square Jt is only contractible on its ambient side
@test pullback(ω1_3d, Jt_emb) isa CoVectorValue{cov2}

# ── Pushforward (value level) ────────────────────────────────────────────────

# The pushforward takes the same Jt = ∇φ as the pullback, but contracts its
# domain (first, covariant) index directly against the K-vector, instead of
# its ambient (second, contravariant) index against a K-form.

e1 = ContraVectorValue{2}(1.0, 0.0)
e2 = ContraVectorValue{2}(0.0, 1.0)

# Identity: pushforward = identity
@test pushforward(e1, Jt_id) == e1
@test pushforward(e2, Jt_id) == e2

biv2 = ExteriorVectorValue{2,2}((1.0,))
@test pushforward(biv2, Jt_id).data == (1.0,)

# Scaling map diag(a,b): pushforward scales each component
@test collect(pushforward(e1, Jt_sc).data) ≈ [a, 0.0]
@test collect(pushforward(e2, Jt_sc).data) ≈ [0.0, b]

# pushforward of ∂1∧∂2 under Jt_sc = diag(a,b): (φ_*biv)^{12} = det(Jt) = a*b
@test collect(pushforward(biv2, Jt_sc).data) ≈ [a*b]   atol=1e-12

# Shear ℝ² → ℝ²: Jt = [1 1; 0 1] is square but not symmetric, so unlike the
# diagonal and identity cases above, this pins which index of Jt contracts
# with the K-vector: (φ_*v)^i = Σ_j Jt[j,i] v^j
Jt_sh = TensorValue{cov2,con2,Float64,4}(1.0,0.0,  1.0,1.0)
@test collect(pushforward(e1, Jt_sh).data) ≈ [1.0, 1.0]
@test collect(pushforward(e2, Jt_sh).data) ≈ [0.0, 1.0]

# Embedding ℝ² → ℝ³: pushforward lifts into ℝ³
pf_e1 = pushforward(e1, Jt_emb)
pf_e2 = pushforward(e2, Jt_emb)
@test typeof(pf_e1) <: ContraVectorValue{con3}
@test collect(pf_e1.data) ≈ [1.0, 0.0, 0.0]
@test collect(pf_e2.data) ≈ [0.0, 1.0, 0.0]

# pushforward of ∂1∧∂2 under the embedding: only the (1,2) ambient component
# survives, since Jt's two rows span exactly the first two ambient directions
pf_biv = pushforward(biv2, Jt_emb)
@test typeof(pf_biv) <: ExteriorVectorValue{2,con3}
@test collect(pf_biv.data) ≈ [1.0, 0.0, 0.0]

# The domain index of Jt is the covariant one, so it is the index that
# contracts with the K-vector. Swapping the two rejects the Jacobian rather
# than silently pushing forward along the transposed map.
@test_throws ArgumentError pushforward(e1, transpose(Jt_emb))
@test_throws ArgumentError pushforward(e1, TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0))
@test_throws DimensionMismatch pushforward(e1, TensorValue{cov3,con2,Float64,6}(1.,0,0,1,0,0))

# ── Top-degree and 2-forms/2-vectors in 3D ───────────────────────────────────
# A fully generic Jt3 (all nine entries distinct, none zero) and forms/vectors
# with distinct nonzero components in every slot, so no permutation term of
# the minor expansions below can accidentally vanish and go untested.

Jt3 = TensorValue{cov3,con3,Float64,9}(1.0,4.0,6.0,  2.0,5.0,8.0,  3.0,7.0,9.0)
# Jt3 = [1 2 3; 4 5 7; 6 8 9], det(Jt3) = 7

vol3    = ExteriorFormValue{3,3}((2.5,))
top3vec = ExteriorVectorValue{3,3}((-1.5,))

# K = D = 3: pullback/pushforward both reduce to a single determinant, the
# only case here that exercises the 3-term (K=3) Leibniz expansion
@test collect(pullback(vol3, Jt3).data)       ≈ [2.5*7]    atol=1e-12
@test collect(pushforward(top3vec, Jt3).data) ≈ [-1.5*7]   atol=1e-12

# K = 2 in D = 3: three independent components on each side, so every entry
# of Jt3 and every one of its 2×2 minors participates in some output component
ω2 = ExteriorFormValue{2,3}((1.0, 2.0, 3.0))    # ω_{12}, ω_{13}, ω_{23}
v2 = ExteriorVectorValue{2,3}((4.0, 5.0, 6.0))  # v^{12}, v^{13}, v^{23}

@test collect(pullback(ω2, Jt3).data)    ≈ [-16.0, -40.0, -43.0]   atol=1e-12
@test collect(pushforward(v2, Jt3).data) ≈ [-20.0, -101.0, -100.0] atol=1e-12

# ── Generic (0,K)-tensor pullback (non-alternating) ─────────────────────────
# `pullback` on a plain covariant MultiValue (not an ExteriorFormValue)
# contracts every index densely: (φ*a)_{j1,…,jK} = Σ a_{i1,…,iK} Jt[j1,i1]⋯Jt[jK,iK],
# with no symmetry of `a` assumed or preserved, unlike the alternating-tensor
# `pullback` above.

# Identity: pullback of a generic (non-symmetric) 2-tensor is unchanged,
# just retyped as a dense TensorValue
A2 = TensorValue{cov2,cov2,Float64,4}(1.0, 3.0, 2.0, 4.0)   # A2 = [1 2; 3 4]
pbA2_id = pullback(A2, Jt_id)
@test typeof(pbA2_id) <: TensorValue{cov2,cov2}
@test collect(pbA2_id.data) ≈ collect(A2.data)   atol=1e-12

# Scaling diag(a,b): (φ*A)_{j1,j2} = A_{j1,j2} · Jt[j1,j1] · Jt[j2,j2]
pbA2_sc = pullback(A2, Jt_sc)
@test collect(pbA2_sc.data) ≈ [1.0*a*a, 3.0*b*a, 2.0*a*b, 4.0*b*b]   atol=1e-12

# Fully generic non-diagonal Jacobian: matches the matrix congruence Jt·A·Jtᵀ
Jt_gen = TensorValue{cov2,con2,Float64,4}(1.0,3.0, 2.0,4.0)   # Jt_gen = [1 2; 3 4]
A2g    = TensorValue{cov2,cov2,Float64,4}(5.0,7.0, 6.0,8.0)   # A2g = [5 6; 7 8]
pbA2g  = pullback(A2g, Jt_gen)
# Jt_gen·A2g·Jt_genᵀ = [63 145; 143 329]
@test collect(pbA2g.data) ≈ [63.0, 143.0, 145.0, 329.0]   atol=1e-12

# Embedding ℝ² → ℝ³: pulling back a 2-tensor is exactly its top-left 2×2 block
A3 = TensorValue{cov3,cov3,Float64,9}(11.,21.,31., 12.,22.,32., 13.,23.,33.)
pbA3 = pullback(A3, Jt_emb)
@test typeof(pbA3) <: TensorValue{cov2,cov2}
@test collect(pbA3.data) ≈ [11.0, 21.0, 12.0, 22.0]   atol=1e-12

# 3rd order tensor: a non-trivial permutation Jacobian only ever selects a
# single term per output component, cross-checking the K ≥ 3
# (NTensorValue) index bookkeeping without a large hand-worked
# contraction. π(j) = perm[j].
perm    = (2, 3, 1)
Jt_perm = TensorValue{cov3,con3,Float64,9}(0.,0.,1.,  1.,0.,0.,  0.,1.,0.)
Araw = [100i + 10j + k for i in 1:3, j in 1:3, k in 1:3]
a3D  = NTensorValue{Tuple{cov3,cov3,cov3}}(Float64.(Araw))
pb3D = pullback(a3D, Jt_perm)
@test typeof(pb3D) <: NTensorValue{Tuple{cov3,cov3,cov3}}
expected3D = [Araw[perm[j1],perm[j2],perm[j3]] for j1 in 1:3, j2 in 1:3, j3 in 1:3]
@test collect(Tuple(pb3D)) ≈ Float64.(vec(expected3D))   atol=1e-12

# Errors mirror the alternating-tensor pullback's guards, plus rejecting an
# un-varianced (plain `Int`-dimensioned) index of `a` (undefined: it is neither
# covariant nor contravariant) and a dimension mismatch against Jt's image
@test_throws ArgumentError pullback(A2, transpose(Jt_emb))
@test_throws ArgumentError pullback(A2, TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0))
@test_throws DimensionMismatch pullback(A2, TensorValue{cov2,con3,Float64,6}(1.,0,0,1,0,0))
@test_throws ArgumentError pullback(TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0), Jt_id)

# ── Mixed variance (r,s,0)-tensor pullback ───────────────────────────────────
# A contravariant index of `a` is contracted with `iJt = pinvJt(Jt)` instead of
# `Jt`; for a square Jt this is exact, for an embedding (Dm < Dn) it is only
# exact for a contravariant index that truly lies in the image tangent space.

# pinvJt itself: for a square varianced Jt, pinvJt(Jt) == inv(Jt), and both
# Jt⋅pinvJt(Jt) and pinvJt(Jt)⋅Jt are the identity
iJt_id = pinvJt(Jt_gen)
@test collect(iJt_id.data) ≈ collect(inv(Jt_gen).data)   atol=1e-12
@test collect((Jt_gen ⋅ iJt_id).data) ≈ [1.0,0.0,0.0,1.0]   atol=1e-12
@test collect((iJt_id ⋅ Jt_gen).data) ≈ [1.0,0.0,0.0,1.0]   atol=1e-12

# for an embedding (Dm < Dn), pinvJt(Jt) should be a left inverse: Jt⋅pinvJt(Jt)
# is the identity on the domain, pinvJt(Jt)⋅Jt a projector on the image.
@test collect((Jt_emb ⋅ pinvJt(Jt_emb)).data) ≈ [1.0,0.0,0.0,1.0]   atol=1e-12

# transpose only reorders the two indices, it does not flip Cov/Con, so
# transpose(Jt_emb) is still a valid (opposite-variance) argument for pinvJt,
# just shaped like the Jacobian of a map the other way round (Dm=3 > Dn=2).
# pinvJt then commutes with transpose exactly, as for inv: (Aᵀ)⁺ = (A⁺)ᵀ.
@test collect(pinvJt(transpose(Jt_emb)).data) ≈ collect(transpose(pinvJt(Jt_emb)).data)   atol=1e-12

# pinvJt returns the Moore-Penrose pseudo-inverse of Jt for any rank. The closed
# forms above are only used when Jt has maximum rank, that is
# rank(Jt) == min(D1,D2), because they are two orders of magnitude faster than
# the singular value decomposition pinvJt falls back to otherwise.
Jt_sq_def   = TensorValue{cov2,con2,Float64,4}(1.0,2.0, 2.0,4.0)           # rank 1, D1==D2
Jt_wide_def = TensorValue{cov2,con3,Float64,6}(1.0,2.0, 2.0,4.0, 3.0,6.0)  # rank 1, D1<D2
Jt_tall_def = TensorValue{cov3,con2,Float64,6}(1.0,2.0,3.0, 2.0,4.0,6.0)   # rank 1, D1>D2

# both paths agree with LinearAlgebra.pinv, the branch is an optimisation only
@test collect(pinvJt(Jt_gen).data)      ≈ vec(pinv([1.0 2.0; 3.0 4.0]))           atol=1e-12
@test collect(pinvJt(Jt_emb).data)      ≈ vec(pinv([1.0 0.0 0.0; 0.0 1.0 0.0]))   atol=1e-12
@test collect(pinvJt(Jt_sq_def).data)   ≈ vec(pinv([1.0 2.0; 2.0 4.0]))           atol=1e-12
@test collect(pinvJt(Jt_wide_def).data) ≈ vec(pinv([1.0 2.0 3.0; 2.0 4.0 6.0]))   atol=1e-12
@test collect(pinvJt(Jt_tall_def).data) ≈ vec(pinv([1.0 2.0; 2.0 4.0; 3.0 6.0]))  atol=1e-12

# the four Moore-Penrose relations, on a rank deficient Jt where neither
# Jt⋅pinvJt(Jt) nor pinvJt(Jt)⋅Jt is an identity
iJt_sq_def = pinvJt(Jt_sq_def)
@test collect((iJt_sq_def ⋅ Jt_sq_def ⋅ iJt_sq_def).data) ≈ collect(iJt_sq_def.data)   atol=1e-12
@test collect((Jt_sq_def ⋅ iJt_sq_def ⋅ Jt_sq_def).data) ≈ collect(Jt_sq_def.data)   atol=1e-12
@test collect((Jt_sq_def ⋅ iJt_sq_def).data) ≈ collect(transpose(Jt_sq_def ⋅ iJt_sq_def).data)   atol=1e-12
@test collect((iJt_sq_def ⋅ Jt_sq_def).data) ≈ collect(transpose(iJt_sq_def ⋅ Jt_sq_def).data)   atol=1e-12

iJt_wide_def = pinvJt(Jt_wide_def)
@test collect((iJt_wide_def ⋅ Jt_wide_def ⋅ iJt_wide_def).data) ≈ collect(iJt_wide_def.data)   atol=1e-12
@test collect((Jt_wide_def ⋅ iJt_wide_def ⋅ Jt_wide_def).data) ≈ collect(Jt_wide_def.data)   atol=1e-12

# transpose commutes with pinvJt for a rank deficient Jt too
@test collect(pinvJt(transpose(Jt_tall_def)).data) ≈ collect(transpose(pinvJt(Jt_tall_def)).data)   atol=1e-12

# the gradient of a constant map is zero, and its pseudo-inverse is zero, where
# the closed form inv(Jt) is singular
@test pinvJt(zero(TensorValue{cov2,con2,Float64,4})) === zero(TensorValue{cov2,con2,Float64,4})

# the fallback preserves the argument eltype, as the closed forms do
@test eltype(pinvJt(TensorValue{cov2,con2,Float32,4}(1.0f0,2.0f0, 2.0f0,4.0f0))) === Float32
@test eltype(pinvJt(TensorValue{cov2,con2,Int,4}(1,2, 2,4))) === Float64

# For a *square* Jt, the metric drops out of the formula algebraically
# (g⋅Jtᵀ⋅inv(Jt⋅g⋅Jtᵀ) = inv(Jt) for any invertible g, since Jt is invertible
# too), so pinvJt's plain `inv(Jt)` already agrees with the metric-aware
# formula for any metric.
Jt_sq_num = [1.0 2.0; 3.0 4.0]
g_sq_num  = [2.0 0.3; 0.3 1.5]
invJt_sq_ref_num = g_sq_num * Jt_sq_num' * inv(Jt_sq_num * g_sq_num * Jt_sq_num')
Jt_sq_ex     = TensorValue{cov2,con2,Float64,4}(vec(Jt_sq_num)...)
invJt_sq_ref = TensorValue{cov2,con2,Float64,4}(vec(invJt_sq_ref_num)...)
@test collect(pinvJt(Jt_sq_ex).data) ≈ collect(invJt_sq_ref.data)   atol=1e-12

# mixed Kronecker delta δ^i_j pulls back to itself, for any invertible Jt
δ = TensorValue{con2,cov2,Float64,4}(1.0,0.0,0.0,1.0)
pbδ = pullback(δ, Jt_gen)
@test typeof(pbδ) <: TensorValue{con2,cov2}
@test collect(pbδ.data) ≈ [1.0,0.0,0.0,1.0]   atol=1e-12

# (1,1) tensor (linear map): pullback is the similarity transform iJtᵀ·A·Jtᵀ
# (independently derived: push a domain vector forward with Jt, apply A in the
# image, then pull the result back with iJt = pinvJt(Jt))
Amix = TensorValue{con2,cov2,Float64,4}(5.0,7.0,6.0,8.0)   # Amix = [5 6; 7 8]
pbAmix = pullback(Amix, Jt_gen)
Jtmat_gen, Amat_mix = [1.0 2.0; 3.0 4.0], [5.0 6.0; 7.0 8.0]
expected_mix = inv(Jtmat_gen)' * Amat_mix * Jtmat_gen'
@test [pbAmix[1,1] pbAmix[1,2]; pbAmix[2,1] pbAmix[2,2]] ≈ expected_mix   atol=1e-12

# embedding case: an ambient contravariant vector that truly lies in the image
# tangent space of `Jt_emb` (ℝ²→ℝ³) is pulled back exactly
w = ContraVectorValue{con2}(2.0,5.0)
v_ambient = ContraVectorValue{con3}(2.0,5.0,0.0)   # = Jt_embᵀ⋅w, in the image tangent space
# pulling back a contravariant index needs pinvJt(Jt_emb) internally
@test typeof(pullback(v_ambient, Jt_emb)) <: ContraVectorValue{con2}
@test collect(Tuple(pullback(v_ambient, Jt_emb))) ≈ collect(Tuple(w))   atol=1e-12

# eltype promotion: an Int `a` pulled back through a contravariant index needs
# the (Int/Int) division inside pinvJt, so the result is not Int
v_int = ContraVectorValue{con2}(1,1)
Jt_int = TensorValue{cov2,con2,Int,4}(2,0,0,2)
@test eltype(pullback(v_int, Jt_int)) == Float64

# a contravariant index with Dm > Dn is not supported (pinvJt has no right inverse)
Jt_wide = TensorValue{cov3,con2,Float64,6}(1.,0.,0., 0.,1.,0.)
@test_throws ArgumentError pullback(w, Jt_wide)

# ── Generic (r,s,0)-tensor pushforward ──────────────────────────────────────
# The exact mirror of the generic pullback: a contravariant index is contracted
# with `Jt`, a covariant one with `iJt = pinvJt(Jt)`, and no symmetry of `a` is
# assumed or preserved.

A2con = TensorValue{con2,con2,Float64,4}(1.0, 3.0, 2.0, 4.0)   # A2con = [1 2; 3 4]
@test typeof(pushforward(A2con, Jt_id)) <: TensorValue{con2,con2}
@test collect(pushforward(A2con, Jt_id).data) ≈ collect(A2con.data)   atol=1e-12

# Scaling diag(a,b): (φ_*A)^{i1,i2} = A^{i1,i2} · Jt[i1,i1] · Jt[i2,i2]
@test collect(pushforward(A2con, Jt_sc).data) ≈ [1.0*a*a, 3.0*b*a, 2.0*a*b, 4.0*b*b] atol=1e-12

# Fully generic non-diagonal Jacobian: a twice contravariant index pair gives
# the congruence Jtᵀ·A·Jt, the transpose of the twice covariant Jt·A·Jtᵀ
A2gcon = TensorValue{con2,con2,Float64,4}(5.0,7.0, 6.0,8.0)   # A2gcon = [5 6; 7 8]
pfA2g  = pushforward(A2gcon, Jt_gen)
A2g_num = [5.0 6.0; 7.0 8.0]
@test [pfA2g[1,1] pfA2g[1,2]; pfA2g[2,1] pfA2g[2,2]] ≈
  Jtmat_gen' * A2g_num * Jtmat_gen   atol=1e-12

# Embedding ℝ² → ℝ³: the pushed forward 2-tensor is the top-left 2×2 block of
# the result, the ambient directions outside the image carrying nothing
A2emb = pushforward(A2gcon, Jt_emb)
@test typeof(A2emb) <: TensorValue{con3,con3}
@test collect(A2emb.data) ≈ [5.0,7.0,0.0, 6.0,8.0,0.0, 0.0,0.0,0.0]   atol=1e-12

# a mixed (1,1) tensor pushes its contravariant index forward with Jt and its
# covariant one with iJt, giving the similarity transform Jtᵀ·A·iJtᵀ
Amixcon = TensorValue{con2,cov2,Float64,4}(5.0,7.0, 6.0,8.0)
pfAmix  = pushforward(Amixcon, Jt_gen)
@test [pfAmix[1,1] pfAmix[1,2]; pfAmix[2,1] pfAmix[2,2]] ≈
  Jtmat_gen' * Amat_mix * inv(Jtmat_gen)'   atol=1e-12

# 3rd order tensor under the permutation Jacobian, the mirror of the pullback
# check above: the pushforward contracts Jt's other index, so it selects the
# inverse permutation. π(j) = perm[j], π⁻¹(i) = iperm[i].
a3Dcon = NTensorValue{Tuple{con3,con3,con3}}(Float64.(Araw))
pf3D   = pushforward(a3Dcon, Jt_perm)
@test typeof(pf3D) <: NTensorValue{Tuple{con3,con3,con3}}
iperm        = (3, 1, 2)
expected3Dpf = [Araw[iperm[i1],iperm[i2],iperm[i3]] for i1 in 1:3, i2 in 1:3, i3 in 1:3]
@test collect(Tuple(pf3D)) ≈ Float64.(vec(expected3Dpf))   atol=1e-12

# Errors mirror the generic pullback's guards
@test_throws ArgumentError pushforward(A2con, transpose(Jt_emb))
@test_throws ArgumentError pushforward(A2con, TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0))
@test_throws DimensionMismatch pushforward(A2con, TensorValue{cov3,con2,Float64,6}(1.,0,0,1,0,0))
@test_throws ArgumentError pushforward(TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0), Jt_id)

# ── Transport against the direction of φ ────────────────────────────────────
# Pushing a K-form forward and pulling a K-vector back are the transports by
# ψ = φ⁻¹, computed from iJt = pinvJt(Jt). All four alternating transports keep
# the alternating type, which the dense (r,s,0) methods above would drop since
# an ExteriorFormValue is a MultiValue too.

iJt_gen = pinvJt(Jt_gen)
@test typeof(pushforward(vol2, Jt_gen))          <: ExteriorFormValue{2,cov2}
@test typeof(pullback(biv2, Jt_gen))             <: ExteriorVectorValue{2,con2}
@test typeof(pushforward(vol2, Jt_gen, iJt_gen)) <: ExteriorFormValue{2,cov2}
@test typeof(pullback(biv2, Jt_gen, iJt_gen))    <: ExteriorVectorValue{2,con2}
@test typeof(pullback(vol2, Jt_gen, iJt_gen))    <: ExteriorFormValue{2,cov2}
@test typeof(pushforward(biv2, Jt_gen, iJt_gen)) <: ExteriorVectorValue{2,con2}

# transporting the volume form against φ scales it by det(iJt) = 1/det(Jt)
@test collect(pushforward(vol2, Jt_gen).data) ≈ [1.0/det(Jtmat_gen)]   atol=1e-12
@test collect(pullback(biv2, Jt_gen).data)    ≈ [1.0/det(Jtmat_gen)]   atol=1e-12

# round trip on a square Jt: the transports invert each other pairwise
@test collect(Tuple(pullback(pushforward(ω1, Jt_gen), Jt_gen))) ≈ collect(Tuple(ω1)) atol=1e-12
@test collect(Tuple(pushforward(pullback(e1, Jt_gen), Jt_gen))) ≈ collect(Tuple(e1)) atol=1e-12
@test collect(pullback(pushforward(vol2, Jt_gen), Jt_gen).data) ≈ collect(vol2.data) atol=1e-12
@test collect(pushforward(pullback(biv2, Jt_gen), Jt_gen).data) ≈ collect(biv2.data) atol=1e-12
@test collect(pullback(pushforward(A2, Jt_gen), Jt_gen).data)   ≈ collect(A2.data)   atol=1e-12
@test collect(pullback(pushforward(Amix, Jt_gen), Jt_gen).data) ≈ collect(Amix.data) atol=1e-12

# the three argument forms only save recomputing iJt = pinvJt(Jt)
@test pushforward(ω1,   Jt_gen, iJt_gen) === pushforward(ω1,   Jt_gen)
@test pushforward(e1,   Jt_gen, iJt_gen) === pushforward(e1,   Jt_gen)
@test pushforward(A2,   Jt_gen, iJt_gen) === pushforward(A2,   Jt_gen)
@test pushforward(Amix, Jt_gen, iJt_gen) === pushforward(Amix, Jt_gen)
@test pullback(ω1,   Jt_gen, iJt_gen) === pullback(ω1,   Jt_gen)
@test pullback(e1,   Jt_gen, iJt_gen) === pullback(e1,   Jt_gen)
@test pullback(A2,   Jt_gen, iJt_gen) === pullback(A2,   Jt_gen)
@test pullback(Amix, Jt_gen, iJt_gen) === pullback(Amix, Jt_gen)

# pinvJt swaps the position and the variance of both indices, so a single shape
# can be pinvJt(Jt) and every three argument form rejects any other
@test_throws ArgumentError pullback(A2,      Jt_gen, transpose(iJt_gen))
@test_throws ArgumentError pushforward(A2con, Jt_gen, transpose(iJt_gen))
@test_throws ArgumentError pullback(vol2,    Jt_gen, transpose(iJt_gen))
@test_throws ArgumentError pushforward(biv2, Jt_gen, transpose(iJt_gen))
@test_throws ArgumentError pushforward(vol2, Jt_gen, transpose(iJt_gen))
@test_throws ArgumentError pullback(biv2,    Jt_gen, transpose(iJt_gen))

# Dm > Dn rejects exactly the transports that need iJt: a contravariant index
# still pushes forward, and a covariant index still pulls back, with Jt alone
@test pushforward(ContraVectorValue{con3}(1.0,2.0,3.0), Jt_wide) isa ContraVectorValue{con2}
@test pushforward(ExteriorVectorValue{2,3}((1.,2.,3.)), Jt_wide) isa ExteriorVectorValue{2,con2}
@test pushforward(TensorValue{con3,con3,Float64,9}(1.,0,0, 0,1.,0, 0,0,1.), Jt_wide) isa TensorValue{con2,con2}
@test pullback(CoVectorValue{cov2}(1.0,2.0), Jt_wide) isa CoVectorValue{cov3}
@test_throws ArgumentError pushforward(CoVectorValue{cov3}(1.0,2.0,3.0), Jt_wide)
@test_throws ArgumentError pushforward(ExteriorFormValue{2,3}((1.,2.,3.)), Jt_wide)
@test_throws ArgumentError pushforward(TensorValue{cov3,cov3,Float64,9}(1.,0,0, 0,1.,0, 0,0,1.), Jt_wide)
@test_throws ArgumentError pullback(ExteriorVectorValue{2,2}((1.0,)), Jt_wide)

# the eltype only promotes against the matrices the indices of `a` reach: `Jt`
# for a contravariant index pushed forward, `iJt` for a covariant one, and both
# for a mixed tensor
Jt_32  = TensorValue{cov2,con2,Float32,4}(1f0,0f0,0f0,1f0)
iJt_64 = TensorValue{cov2,con2,Float64,4}(1.0,0.0,0.0,1.0)
Amix_i = TensorValue{con2,cov2,Int,4}(1,0,0,1)
@test eltype(pushforward(ContraVectorValue{con2,Int}((1,2)), Jt_32, iJt_64)) === Float32
@test eltype(pushforward(CoVectorValue{cov2,Int}((1,2)),     Jt_32, iJt_64)) === Float64
@test eltype(pushforward(Amix_i, Jt_32, iJt_64))                             === Float64
@test eltype(pullback(CoVectorValue{cov2,Int}((1,2)),     Jt_32, iJt_64))    === Float32
@test eltype(pullback(ContraVectorValue{con2,Int}((1,2)), Jt_32, iJt_64))    === Float64
@test eltype(pullback(Amix_i, Jt_32, iJt_64))                                === Float64

# ── A ProductTensorValue jacobian  ─────────────────────────────────

Jt_pt  = ProductTensorValue{SymNTensorValue{Tuple{cov2}},ContraVectorValue{con2}}(1.0,3.0,2.0,4.0)
iJt_pt = ProductTensorValue{SymNTensorValue{Tuple{cov2}},ContraVectorValue{con2}}(Tuple(iJt_gen))
@test Jt_pt isa MultiValue{Tuple{cov2,con2}}
@test get_array(Jt_pt) == get_array(Jt_gen)

for x in (ω1, e1, vol2, biv2, A2, A2con, Amix)
  @test pullback(x, Jt_pt)             === pullback(x, Jt_gen)
  @test pushforward(x, Jt_pt)          === pushforward(x, Jt_gen)
  @test pullback(x, Jt_pt, iJt_pt)     === pullback(x, Jt_gen, iJt_gen)
  @test pushforward(x, Jt_pt, iJt_pt)  === pushforward(x, Jt_gen, iJt_gen)
end

# A dense jacobian pairs with a factored pseudo-inverse, and conversely
@test pullback(A2, Jt_gen, iJt_pt)        === pullback(A2, Jt_gen, iJt_gen)
@test pushforward(A2con, Jt_pt, iJt_gen)  === pushforward(A2con, Jt_gen, iJt_gen)

# A product tensor transported by a product tensor keeps its factors, as the
# carrier of the jacobian does not reach the transport of the value
∂ω = product_tensor(SymNTensorValue{Tuple{cov2}}(1.0,2.0), CoVectorValue{cov2}(3.0,4.0))
@test pullback(∂ω, Jt_pt) isa
      ProductTensorValue{SymNTensorValue{Tuple{cov2}},CoVectorValue{cov2}}
@test pullback(∂ω, Jt_pt)            === pullback(∂ω, Jt_gen)
@test pushforward(∂ω, Jt_pt)         === pushforward(∂ω, Jt_gen)
@test pullback(∂ω, Jt_pt, iJt_pt)    === pullback(∂ω, Jt_gen, iJt_gen)
@test pushforward(∂ω, Jt_pt, iJt_pt) === pushforward(∂ω, Jt_gen, iJt_gen)

# ── Koszul operator κ_x(ω) = ι_x(ω) ─────────────────────────────────────────

# κ on a 1-form: at x₀ = (2,3), ω = x¹dx¹ + x²dx² evaluates to (2,3);
# κ_x₀(ω) = x¹·x¹ + x²·x² = 4 + 9 = 13
x₀  = VectorValue(2.0, 3.0)
ωx₀ = CoVectorValue{2}((2.0, 3.0))
κω1 = koszul(x₀, ωx₀)
@test typeof(κω1) <: ExteriorFormValue{0,cov2}
@test κω1.data[1] ≈ 13.0   atol=1e-12

# κ on a 2-form: ι_{(2,3)}(x¹ dx¹∧dx²) at x₀=(2,3), coefficient x¹=2:
# d[1] = -x₀[2]·2 = -6, d[2] = x₀[1]·2 = 4
ω2x₀ = ExteriorFormValue{2,2}((2.0,))
κω2  = koszul(x₀, ω2x₀)
@test typeof(κω2) <: CoVectorValue{cov2}
@test collect(κω2.data) ≈ [-6.0, 4.0]   atol=1e-12

# koszul == interior_product (numeric, ambient/barycentric indices)
ω3 = CoVectorValue{3}((0.5, 0.3, 0.2))
v3 = VectorValue(0.5, 0.3, 0.2)
@test koszul(v3, ω3).data[1] ≈ interior_product(v3, ω3).data[1] ≈ 0.38

# ι on an ambient 2-form: ι_{(0.5,0.3,0.2)}(dλ¹∧dλ²) = 0.5 dλ² − 0.3 dλ¹
e12 = ExteriorFormValue{2,3}((1.0, 0.0, 0.0))
k = koszul(VectorValue(0.5, 0.3, 0.2), e12)
@test typeof(k) <: CoVectorValue{cov3}
@test collect(k.data) ≈ [-0.3, 0.5, 0.0]   atol=1e-12

end # module
