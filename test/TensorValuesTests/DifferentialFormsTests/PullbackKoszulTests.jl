module PullbackKoszulTests
# Value-level pullback/pushforward under a linear map (Jacobian), the Koszul
# operator κ_x(ω) = ι_x(ω).

using Gridap.TensorValues
using Test

# ── Pullback (value level) ───────────────────────────────────────────────────

# Identity map on ℝ²: J = I₂, pullback should be identity
J_id = TensorValue{2,2,Float64,4}(1.0,0.0,0.0,1.0)
ω1   = DifferentialFormValue{1,2}((3.0, 5.0))
@test pullback(ω1, J_id).data == (3.0, 5.0)

vol2 = DifferentialFormValue{2,2}((1.0,))
@test pullback(vol2, J_id).data == (1.0,)

# Scaling map ℝ² → ℝ²: φ(ξ) = (aξ¹, bξ²), J = diag(a,b)
a, b = 2.0, 3.0
J_sc = TensorValue{2,2,Float64,4}(a,0.0,0.0,b)

# pullback of dx¹: (φ*dx¹)_j = J[1,j] → a dx¹  (only j=1 gives non-zero)
@test collect(pullback(DifferentialFormValue{1,2}((1.0,0.0)), J_sc).data) ≈ [a, 0.0]
@test collect(pullback(DifferentialFormValue{1,2}((0.0,1.0)), J_sc).data) ≈ [0.0, b]

# pullback of dx¹∧dx² under J_sc = diag(a,b): (φ*vol)_12 = det(J) = a*b
@test collect(pullback(vol2, J_sc).data) ≈ [a*b]   atol=1e-12

# Embedding ℝ² → ℝ³: project onto (x¹,x²,0) plane
# J = [[1,0],[0,1],[0,0]] (3 rows, 2 cols)
J_emb = TensorValue{3,2,Float64,6}(1.0,0.0,0.0,  0.0,1.0,0.0)
ω1_3d = DifferentialFormValue{1,3}((1.0, 0.0, 0.0))  # dx¹ in ℝ³
pb    = pullback(ω1_3d, J_emb)
@test typeof(pb) <: DifferentialFormValue{1,2}
@test pb.data == (1.0, 0.0)   # dx¹ pulls back to dx¹ in ℝ²

ω1_3d_2 = DifferentialFormValue{1,3}((0.0, 1.0, 0.0))  # dx² in ℝ³
@test pullback(ω1_3d_2, J_emb).data == (0.0, 1.0)       # dx² in ℝ²

ω1_3d_3 = DifferentialFormValue{1,3}((0.0, 0.0, 1.0))  # dx³ in ℝ³ → 0 in ℝ²
@test collect(pullback(ω1_3d_3, J_emb).data) ≈ [0.0, 0.0]

# ── Pushforward (value level) ────────────────────────────────────────────────

e1 = VectorValue(1.0, 0.0)
e2 = VectorValue(0.0, 1.0)

# Identity: pushforward = identity
@test pushforward(e1, J_id) == e1
@test pushforward(e2, J_id) == e2

# Scaling map diag(a,b): pushforward scales each component
@test collect(pushforward(e1, J_sc).data) ≈ [a, 0.0]
@test collect(pushforward(e2, J_sc).data) ≈ [0.0, b]

# Embedding ℝ² → ℝ³: pushforward lifts into ℝ³
pf_e1 = pushforward(e1, J_emb)
@test typeof(pf_e1) <: VectorValue{3}
@test collect(pf_e1.data) ≈ [1.0, 0.0, 0.0]

# ── Koszul operator κ_x(ω) = ι_x(ω) ─────────────────────────────────────────

# κ on a 1-form: at x₀ = (2,3), ω = x¹dx¹ + x²dx² evaluates to (2,3);
# κ_x₀(ω) = x¹·x¹ + x²·x² = 4 + 9 = 13
x₀  = VectorValue(2.0, 3.0)
ωx₀ = DifferentialFormValue{1,2}((2.0, 3.0))
κω1 = koszul(x₀, ωx₀)
@test typeof(κω1) <: DifferentialFormValue{0,2}
@test κω1.data[1] ≈ 13.0   atol=1e-12

# κ on a 2-form: ι_{(2,3)}(x¹ dx¹∧dx²) at x₀=(2,3), coefficient x¹=2:
# d[1] = -x₀[2]·2 = -6, d[2] = x₀[1]·2 = 4
ω2x₀ = DifferentialFormValue{2,2}((2.0,))
κω2  = koszul(x₀, ω2x₀)
@test typeof(κω2) <: DifferentialFormValue{1,2}
@test collect(κω2.data) ≈ [-6.0, 4.0]   atol=1e-12

# koszul == interior_product (numeric, ambient/barycentric indices)
ω3 = DifferentialFormValue{1,3}((0.5, 0.3, 0.2))
v3 = VectorValue(0.5, 0.3, 0.2)
@test koszul(v3, ω3).data[1] ≈ interior_product(v3, ω3).data[1] ≈ 0.38

# ι on an ambient 2-form: ι_{(0.5,0.3,0.2)}(dλ¹∧dλ²) = 0.5 dλ² − 0.3 dλ¹
e12 = DifferentialFormValue{2,3}((1.0, 0.0, 0.0))
k = koszul(VectorValue(0.5, 0.3, 0.2), e12)
@test typeof(k) <: DifferentialFormValue{1,3}
@test collect(k.data) ≈ [-0.3, 0.5, 0.0]   atol=1e-12

end # module
