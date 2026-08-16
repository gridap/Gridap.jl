module FormBridgeTests
# Value-level bridges between the form algebra and Gridap's FE machinery:
#   vol_coeff       — scalar content of a top D-form (makes ∫(ω ∧ ⋆η) integrable)
#   to_1form/from_1form, to_Kform/from_Kform, to_0form/to_Dform — VectorValue proxies
#   jac_to_2form    — exterior derivative of a 1-form from its (transposed) Jacobian

using Gridap.TensorValues
using Test

# ── vol_coeff ────────────────────────────────────────────────────────────────

# 2D: single-component 2-form
ω2d = DifferentialFormValue{2,2}((5.0,))
@test vol_coeff(ω2d) === 5.0
@test vol_coeff(ω2d) isa Float64

# 3D: single-component 3-form
ω3d = DifferentialFormValue{3,3}((-2.5,))
@test vol_coeff(ω3d) === -2.5

# vol_coeff extracts the scalar content of ω ∧ ⋆η
# In 2D: ω = (a, b), ⋆ω = (-b, a); ω ∧ ⋆ω = (a²+b²) dx¹∧dx²
ω = DifferentialFormValue{1,2}((3.0, 4.0))
@test vol_coeff(ω ∧ hodge_star(ω)) ≈ 25.0   atol=1e-12  # 3²+4² = 25

# vol_coeff(du ∧ ⋆dv) = ∇u·∇v
du = DifferentialFormValue{1,2}((1.0, 2.0))
dv = DifferentialFormValue{1,2}((3.0, 4.0))
@test vol_coeff(du ∧ hodge_star(dv)) ≈ 11.0   # 1·3 + 2·4

# ── to_1form / from_1form ────────────────────────────────────────────────────

v = VectorValue(1.0, 2.0)
ω1 = to_1form(v)
@test typeof(ω1) <: DifferentialFormValue{1,2}
@test ω1.data == (1.0, 2.0)

v2 = from_1form(ω1)
@test typeof(v2) <: VectorValue{2}
@test v2 == v

# Roundtrip
w = VectorValue(-1.0, 3.0, 0.5)
@test from_1form(to_1form(w)) == w

# 3D
v3 = VectorValue(1.0, 2.0, 3.0)
ω3 = to_1form(v3)
@test typeof(ω3) <: DifferentialFormValue{1,3}
@test from_1form(ω3) == v3

# ── to_Kform / from_Kform: general K ─────────────────────────────────────────

# K=2 in D=3: binomial(3,2) = 3 components
v_2form = VectorValue(1.0, -2.0, 3.0)
ω_2form = to_Kform(v_2form, Val(2), Val(3))
@test typeof(ω_2form) <: DifferentialFormValue{2,3}
@test ω_2form.data == (1.0, -2.0, 3.0)
@test from_Kform(ω_2form) == v_2form

# ── to_0form / to_Dform ──────────────────────────────────────────────────────

@test to_0form(3.0, Val(2)) isa DifferentialFormValue{0,2}
@test to_0form(3.0, Val(2)).data == (3.0,)
@test to_Dform(5.0, Val(2)) isa DifferentialFormValue{2,2}
@test to_Dform(5.0, Val(2)).data == (5.0,)
@test to_Dform(7.0, Val(3)) isa DifferentialFormValue{3,3}

# ── jac_to_2form ─────────────────────────────────────────────────────────────
# Gridap's ∇(u) uses the TRANSPOSED Jacobian: J[i,j] = ∂u_j/∂x_i
# (dω)_{a<b} = ∂ω_b/∂x^a − ∂ω_a/∂x^b = J[a,b] − J[b,a]
# TensorValue{2,2,T,4}(a,b,c,d) stores in column-major:
#   T[1,1]=a, T[2,1]=b, T[1,2]=c, T[2,2]=d

# Identity Jacobian: ω = x¹dx¹+x²dx² has dω = 0 (J = I, symmetric)
J_id = TensorValue{2,2,Float64,4}(1.0, 0.0, 0.0, 1.0)
@test vol_coeff(jac_to_2form(J_id)) == 0.0

# Rotation u = (−x², x¹): Gridap's ∇ gives J[i,j]=∂u_j/∂x_i
# J[2,1]=∂u₁/∂x₂=−1,  J[1,2]=∂u₂/∂x₁=1  → TensorValue(0,−1,1,0)
# (dω)₁₂ = J[1,2]−J[2,1] = 1−(−1) = 2
J_rot = TensorValue{2,2,Float64,4}(0.0, -1.0, 1.0, 0.0)
@test vol_coeff(jac_to_2form(J_rot)) ≈ 2.0

# Shear u = (x², 0): J[2,1]=∂u₁/∂x₂=1, others 0 → TensorValue(0,1,0,0)
# (dω)₁₂ = J[1,2]−J[2,1] = 0−1 = −1
J_sh = TensorValue{2,2,Float64,4}(0.0, 1.0, 0.0, 0.0)
@test vol_coeff(jac_to_2form(J_sh)) ≈ -1.0

# 3D: ω = (x¹, x², x³), J = I₃ → dω = 0 (symmetric, convention-independent)
J3 = TensorValue{3,3,Float64,9}(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0)
dω3 = jac_to_2form(J3)
@test typeof(dω3) <: DifferentialFormValue{2,3}
@test all(iszero, dω3.data)

end # module
